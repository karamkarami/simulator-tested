/*
 * ppi_qt.cpp - interactive Qt PPI (plan position indicator) of the radar receiver.
 *
 * The worker thread posts antenna positions and plots with ppi_post_antenna()
 * and ppi_post_plot(). These calls only store the data under a mutex; a timer
 * in the GUI thread (25 Hz) takes the new data and repaints the widget.
 *
 * Interaction:
 *   mouse wheel            zoom in/out around the cursor (x1 ... x500)
 *   left drag              pan
 *   double click           reset the view (whole instrumented range)
 *   hover over a plot      data tip: range, azimuth, elevation, radial velocity,
 *                          SNR, time, sub-pulse, hits, beam/cell/bin, pulse number
 *   click on a plot        pin its data tip on the display (click again: remove)
 *   right click            options: what is shown (plots / CFAR detections /
 *                          both), short/long-pulse data, range rings, plot
 *                          labels, antenna trail, which plots are shown (until
 *                          the antenna passes again / 2 scans / all), minimum
 *                          SNR, clear the display, remove data tips, reset
 *                          view, pause, save a screenshot
 *   keys                   + - zoom, arrows pan, R/Home reset, C remove the
 *                          data tips, D next display mode, Space pause/resume,
 *                          S save a screenshot
 * The status line shows the range, azimuth and east/north of the cursor.
 *
 * Panel on the right of the window (PpiWindow): display mode (plots after the
 * plot extractor, CFAR detections before it, or both), a button that saves a
 * screenshot of the PPI, and the CFAR threshold [dB] of the running
 * processing: Apply sends the value to the worker (ppi_take_cfar_request),
 * which reports the threshold it uses back (ppi_post_cfar_threshold); Default
 * restores the design value of the configuration.
 *
 * CFAR detections (ppi_post_detections) are the pruned detections of every
 * CPI update, before the plot extractor. Each is drawn as an arc over the
 * azimuth step between CPI updates, so a target becomes an arc about one
 * beamwidth long (one arc per range cell and Doppler bin that is detected).
 *
 * Plot persistence (default): a plot stays until the antenna sweeps over its
 * azimuth again in a later scan, like on a real PPI. For this the antenna
 * rotation is accumulated ("sweep", degrees turned since the first update), so
 * the scan boundary (the start azimuth of the antenna) does not matter.
 *
 * The file does not use Q_OBJECT, so it builds with or without moc, with
 * Qt 5 (>= 5.12) and Qt 6.
 */
#include "radar_rx.h"

#include <QAction>
#include <QActionGroup>
#include <QApplication>
#include <QCloseEvent>
#include <QColor>
#include <QConicalGradient>
#include <QCheckBox>
#include <QComboBox>
#include <QContextMenuEvent>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QDateTime>
#include <QHBoxLayout>
#include <QLabel>
#include <QPainterPath>
#include <QPushButton>
#include <QVBoxLayout>
#include <QFont>
#include <QFontDatabase>
#include <QFontMetrics>
#include <QKeyEvent>
#include <QMenu>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>
#include <QPen>
#include <QPixmap>
#include <QPointF>
#include <QRectF>
#include <QString>
#include <QStringList>
#include <QTimer>
#include <QToolTip>
#include <QWheelEvent>
#include <QWidget>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <deque>
#include <functional>
#include <mutex>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr double kPi          = 3.14159265358979323846;
constexpr double kDeg         = kPi / 180.0;
constexpr size_t kMaxPlots    = 4000;     /* plots kept on the display              */
constexpr size_t kMaxPending  = 20000;    /* plots waiting for the GUI thread       */
constexpr size_t kMaxDets     = 60000;    /* CFAR detections kept on the display    */
constexpr size_t kMaxPendingDets = 60000; /* detections waiting for the GUI thread  */
constexpr int    kRefreshMs   = 40;       /* GUI refresh period (25 Hz)             */
constexpr double kHeaderPx    = 48.0;     /* status text above the display          */
constexpr double kLabelPx     = 28.0;     /* room for the azimuth labels            */
constexpr double kPlotRadius  = 3.5;      /* plot symbol radius [px]                */
constexpr double kPickPx      = 8.0;      /* data tip: max. distance to a plot [px] */
constexpr double kTrailDeg    = 20.0;     /* fading trail behind the antenna line   */
constexpr double kMaxZoom     = 500.0;

const QColor kColorShort(0, 230, 255);    /* cyan                                   */
const QColor kColorLong(255, 230, 0);     /* yellow                                 */
const QColor kColorDetShort(0, 150, 180); /* CFAR detections: darker cyan / yellow  */
const QColor kColorDetLong(200, 150, 0);
const QColor kColorAntenna(0, 220, 0);    /* green                                  */
const QColor kColorGrid(70, 70, 70);
const QColor kColorRim(120, 120, 120);
const QColor kColorLabel(150, 150, 150);
const QColor kColorPin(80, 230, 80);
const QColor kColorTrack(255, 90, 255);   /* confirmed track: magenta triangle      */
const QColor kColorTentative(140, 140, 140);
constexpr double kTrackPx     = 11.0;     /* half size of the track triangle [px]   */
constexpr double kLeader_s    = 30.0;     /* velocity leader of a track [s]         */

double wrap360(double a)
{
    a = std::fmod(a, 360.0);
    return a < 0.0 ? a + 360.0 : a;
}

double wrap180(double a)
{
    a = wrap360(a);
    return a >= 180.0 ? a - 360.0 : a;
}

/* ------------------------------------------------------------------------ */
/* Data shared between the worker thread and the GUI thread                  */
/* ------------------------------------------------------------------------ */

struct AntennaState {
    bool     valid = false;
    double   az_deg = 0.0;
    double   time_s = 0.0;
    uint32_t scan = 0;
    double   sweep_deg = 0.0;   /* rotation accumulated since the first update */
};

struct PlotMark {
    double    east_m = 0.0;     /* ground position                              */
    double    north_m = 0.0;
    int       pulse = PULSE_SHORT;
    double    sweep_deg = 0.0;  /* antenna sweep at which it pointed at the plot */
    uint32_t  scan = 0;         /* antenna scan when the plot arrived           */
    uint64_t  id = 0;           /* unique number (pinned data tips)             */
    bool      isDet = false;    /* CFAR detection (before the extractor), not a plot */
    double    arc_deg = 0.0;    /* detections: azimuth width of the arc         */
    Detection det;              /* the plot or detection itself (data tip)      */
};

struct CfarState {
    bool     known = false;     /* the worker reported its threshold            */
    double   dB = 0.0, pfa = NAN, default_dB = NAN;
    uint64_t serial = 0;        /* counts the reports                           */
    bool     requested = false; /* Apply pressed, not yet taken by the worker   */
    double   request_dB = 0.0;
};

struct Shared {
    std::mutex            mutex;
    AntennaState          antenna;
    std::vector<PlotMark> newPlots;   /* plots not yet taken by the GUI thread */
    std::vector<PlotMark> newDets;    /* CFAR detections not yet taken         */
    bool                  clear = false;   /* a new run started: clear the display */
    uint64_t              nextId = 1;
    CfarState             cfar;
    std::map<uint32_t, Track> trackUpdates;  /* newest report of every track (deleted ones flagged) */
    ProcessorTiming       timing;          /* accumulated block times               */
    std::string           detName = "CFAR", trackName = "Tracking";
    uint64_t              timingSerial = 0;
};

/* Monotonic clock of the GUI thread's display time [s]. */
double gui_seconds()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

Shared &shared()
{
    static Shared s;
    return s;
}

void reset_shared()
{
    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.antenna = AntennaState();
    s.newPlots.clear();
    s.newDets.clear();
    s.clear = true;
    s.cfar = CfarState();
    s.trackUpdates.clear();
    s.timing = ProcessorTiming();
    s.timingSerial = 0;
}

/* Ground position and antenna sweep of a plot or detection (call with the mutex held). */
PlotMark make_mark(const Detection &d, const AntennaState &a, uint64_t id)
{
    const double az = wrap360(d.az_deg);
    double ground_m = d.range_m;
    if (std::isfinite(d.el_deg)) ground_m *= std::cos(d.el_deg * kDeg);
    PlotMark m;
    m.east_m = ground_m * std::sin(az * kDeg);
    m.north_m = ground_m * std::cos(az * kDeg);
    m.pulse = (d.pulse == PULSE_LONG) ? PULSE_LONG : PULSE_SHORT;
    m.det = d;
    /* the azimuth is close to (just behind) the antenna */
    m.sweep_deg = a.sweep_deg + (a.valid ? wrap180(az - a.az_deg) : 0.0);
    m.scan = a.scan;
    m.id = id;
    return m;
}

/* Rotation from the previous antenna update to az_deg. The antenna turns
 * clockwise, so a gap in the data (e.g. sector transmission) is a step forward;
 * a tiny step backwards (jitter) stays a step backwards. When the scan number
 * jumped by more than one, the missing whole turns are added. */
double antenna_step(const AntennaState &prev, double az_deg, uint32_t scan)
{
    double step = wrap360(az_deg - prev.az_deg);
    if (step > 350.0) step -= 360.0;
    const long long jump = (long long)scan - (long long)prev.scan;
    if (jump > 1) step += 360.0 * (double)(jump - 1);
    return step;
}

/* Ring spacing for the visible range: about 3 to 6 rings. */
double ring_step_m(double range_m)
{
    const double km = range_m / 1e3;
    if (km > 500.0) return 100e3;
    if (km > 200.0) return 50e3;
    if (km > 60.0)  return 20e3;
    if (km > 30.0)  return 10e3;
    if (km > 15.0)  return 5e3;
    if (km > 6.0)   return 2e3;
    if (km > 3.0)   return 1e3;
    if (km > 1.5)   return 500.0;
    if (km > 0.6)   return 200.0;
    return 100.0;
}

QString km_text(double m)
{
    return QString::number(m / 1e3, 'g', 4) + QStringLiteral(" km");
}

QString deg_text(double deg, int decimals)
{
    return QString::number(deg, 'f', decimals) + QChar(0x00B0);
}

QFont pixel_font(const QFont &base, int pixels, bool bold = false)
{
    QFont f(base);
    f.setPixelSize(pixels);
    f.setBold(bold);
    return f;
}

/* Data tip text of a plot or of a CFAR detection. */
QString plot_text(const Detection &d, bool isDet = false)
{
    const QString el = std::isfinite(d.el_deg) ? deg_text(d.el_deg, 2) : QStringLiteral("- (no DBF)");
    QStringList rows;
    rows << QStringLiteral("<b>%1, %2 pulse</b>").arg(isDet ? "CFAR detection" : "Plot")
                                                .arg(d.pulse == PULSE_LONG ? "long" : "short")
         << QStringLiteral("Range: %1 km").arg(d.range_m / 1e3, 0, 'f', 3)
         << QStringLiteral("Azimuth: %1").arg(deg_text(d.az_deg, 2))
         << QStringLiteral("Elevation: %1").arg(el)
         << QStringLiteral("Radial velocity: %1 m/s").arg(d.vr_mps, 0, 'f', 2)
         << QStringLiteral("SNR: %1 dB").arg(d.snr_dB, 0, 'f', 1)
         << QStringLiteral("Time: %1 s").arg(d.time_s, 0, 'f', 3);
    if (!isDet) rows << QStringLiteral("CPI hits: %1").arg(d.hits);
    rows << QStringLiteral("Beam %1, range cell %2, Doppler bin %3").arg(d.beam).arg(d.rIdx).arg(d.dIdx)
         << QStringLiteral("Pulse number %1").arg(d.pulseNumber);
    return rows.join(QStringLiteral("<br>"));
}

/* Event positions (the accessors differ between Qt 5 and Qt 6). */
QPointF event_pos(const QMouseEvent *e)
{
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
    return e->position();
#else
    return e->localPos();
#endif
}

QPointF event_pos(const QWheelEvent *e)
{
#if QT_VERSION >= QT_VERSION_CHECK(5, 14, 0)
    return e->position();
#else
    return e->posF();
#endif
}

/* ------------------------------------------------------------------------ */
/* PPI widget                                                                */
/* ------------------------------------------------------------------------ */

enum class History { UntilSwept, TwoScans, All };
enum class DisplayMode { Plots, Detections, Both };   /* after / before the plot extractor */

class PpiWidget : public QWidget {
public:
    PpiWidget(double maxRange_m, volatile sig_atomic_t *stopFlag)
        : maxRange_m_(maxRange_m > 0.0 ? maxRange_m : 60e3), stopFlag_(stopFlag)
    {
        setAttribute(Qt::WA_OpaquePaintEvent);
        setMouseTracking(true);                 /* hover data tips without a button */
        setFocusPolicy(Qt::StrongFocus);        /* keys */
        setMinimumSize(300, 300);
    }

    void refresh();                     /* takes the posted data (GUI thread) */
    /* GUI thread time spent on the display: taking the data and drawing [ms] */
    double guiMs() const { return guiMs_; }
    void addGuiMs(double ms) { guiMs_ += ms; }
    void setFinished(int workerCode)
    {
        finished_ = true;
        workerCode_ = workerCode;
        update();
    }

    /* view and options, also used by the tests */
    void   zoomAt(double factor, QPointF screen);
    void   panPixels(QPointF delta);
    void   resetView();
    double zoom() const { return zoom_; }
    QPointF viewCentre_m() const { return QPointF(viewE_, viewN_); }
    QPointF screenOf(double east_m, double north_m) { layoutDisplay(); return toScreen(east_m, north_m); }
    /* plot or detection under the cursor (nearest within kPickPx), nullptr if none */
    const PlotMark *pick(QPointF screen) const;
    /* track under the cursor (within kPickPx; it wins over a plot that is not 3 px closer) */
    const Track *pickTrack(QPointF screen, double *dist = nullptr) const;
    QString tipText(QPointF screen) const
    {
        uint64_t key;
        return hoverText(screen, &key);
    }
    QString hoverText(QPointF screen, uint64_t *key) const;
    size_t trackCount() const { return tracks_.size(); }
    size_t confirmedCount() const;
    void   setShowTracks(bool on) { showTracks_ = on; update(); }
    bool   showTracks() const { return showTracks_; }
    void   clickAt(QPointF screen);     /* pin / unpin the data tip of a plot or detection */
    size_t pinnedCount() const { return pinned_.size(); }
    size_t plotCount() const { return plots_.size(); }
    size_t detectionCount() const { return dets_.size(); }
    size_t visibleCount() const;        /* plots shown */
    size_t visibleDetectionCount() const;
    void   setDisplayMode(DisplayMode m) { mode_ = m; update(); }
    DisplayMode displayMode() const { return mode_; }
    QString cursorText() const { return cursorText_; }
    void   setShowShort(bool on) { showShort_ = on; update(); }
    void   setShowLong(bool on) { showLong_ = on; update(); }
    void   setMinSnr(double dB) { minSnr_dB_ = dB; update(); }
    void   setHistory(History h) { history_ = h; refresh(); }
    void   setPaused(bool on) { paused_ = on; update(); }
    bool   paused() const { return paused_; }
    void   clearPlots() { plots_.clear(); dets_.clear(); tracks_.clear(); pinned_.clear(); update(); }
    bool   saveSnapshot(const QString &path);
    QString saveSnapshotDialog();       /* asks for a file name; returns it ('' if cancelled) */

protected:
    void paintEvent(QPaintEvent *event) override;
    void closeEvent(QCloseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void mouseDoubleClickEvent(QMouseEvent *event) override;
    void leaveEvent(QEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
    void contextMenuEvent(QContextMenuEvent *event) override;

private:
    double guiMs_ = 0.0;
    void   layoutDisplay();
    void   drawBackground();
    void   drawAntenna(QPainter &p) const;
    void   drawPlots(QPainter &p) const;
    void   drawDetections(QPainter &p) const;
    void   drawTracks(QPainter &p) const;
    void   drawPins(QPainter &p) const;
    void   drawStatus(QPainter &p) const;
    void   drawLegend(QPainter &p) const;
    bool   expired(const PlotMark &m) const;
    bool   shown(const PlotMark &m) const;
    void   updateCursor(QPointF screen);
    void   setView(double zoom, double e, double n);
    double turned(const PlotMark &m) const { return std::fabs(antenna_.sweep_deg - m.sweep_deg); }
    double scale() const { return radius_px_ * zoom_ / maxRange_m_; }   /* pixels per metre */
    QPointF toScreen(double east_m, double north_m) const
    {
        const double k = scale();
        return QPointF(centre_.x() + k * (east_m - viewE_), centre_.y() - k * (north_m - viewN_));
    }
    QPointF toGround(QPointF s) const
    {
        const double k = scale();
        return QPointF(viewE_ + (s.x() - centre_.x()) / k, viewN_ - (s.y() - centre_.y()) / k);
    }
    QPointF polar(double az_deg, double r_m) const
    {
        return toScreen(r_m * std::sin(az_deg * kDeg), r_m * std::cos(az_deg * kDeg));
    }

    double                maxRange_m_;
    volatile sig_atomic_t *stopFlag_;
    AntennaState          antenna_;
    std::deque<PlotMark>  plots_;
    std::deque<PlotMark>  dets_;        /* CFAR detections before the extractor */
    std::map<uint32_t, Track> tracks_;  /* current tracks by number            */
    std::map<uint32_t, Track> incomingTracks_;
    std::vector<PlotMark> incoming_;    /* buffer swapped with Shared::newPlots */
    std::vector<PlotMark> incomingDets_;
    bool                  finished_ = false;
    int                   workerCode_ = 0;
    QPointF               centre_;      /* screen point of the view centre       */
    double                radius_px_ = 1.0;   /* screen radius of maxRange at zoom 1 */
    QRectF                viewRect_;    /* drawing area below the status lines   */
    QPixmap               background_;  /* rings, spokes and labels            */
    bool                  backgroundDirty_ = true;
    QSize                 backgroundSize_;

    /* view */
    double                zoom_ = 1.0;
    double                viewE_ = 0.0, viewN_ = 0.0;   /* ground point at the view centre */
    /* mouse */
    bool                  dragging_ = false, dragMoved_ = false;
    QPointF               dragStart_;
    double                dragE_ = 0.0, dragN_ = 0.0;
    uint64_t              hoverId_ = 0;
    QString               cursorText_;
    std::set<uint64_t>    pinned_;
    /* options */
    bool                  showShort_ = true, showLong_ = true, showRings_ = true;
    bool                  showLabels_ = false, showTrail_ = true, paused_ = false, showTracks_ = true;
    double                minSnr_dB_ = -1e9;
    History               history_ = History::UntilSwept;
    DisplayMode           mode_ = DisplayMode::Plots;
};

void PpiWidget::refresh()
{
    if (paused_) return;                /* the data waits in Shared */
    bool clear = false;
    incoming_.clear();
    incomingDets_.clear();
    {
        Shared &s = shared();
        std::lock_guard<std::mutex> lock(s.mutex);
        incoming_.swap(s.newPlots);     /* the worker keeps the old capacity */
        incomingDets_.swap(s.newDets);
        incomingTracks_.swap(s.trackUpdates);
        s.trackUpdates.clear();
        antenna_ = s.antenna;
        clear = s.clear;
        s.clear = false;
    }
    if (clear) { plots_.clear(); dets_.clear(); tracks_.clear(); pinned_.clear(); }
    for (const auto &u : incomingTracks_) {             /* newest state of every track */
        if (u.second.deleted) tracks_.erase(u.first);
        else tracks_[u.first] = u.second;
    }
    incomingTracks_.clear();
    auto keep = [this](std::deque<PlotMark> &list, const std::vector<PlotMark> &add, size_t maxSize) {
        list.insert(list.end(), add.begin(), add.end());
        list.erase(std::remove_if(list.begin(), list.end(), [this](const PlotMark &m) { return expired(m); }),
                   list.end());
        while (list.size() > maxSize) list.pop_front();
    };
    keep(plots_, incoming_, kMaxPlots);
    keep(dets_, incomingDets_, kMaxDets);
    for (auto it = pinned_.begin(); it != pinned_.end();) {   /* tips of removed plots and detections */
        const uint64_t id = *it;
        auto same = [id](const PlotMark &m) { return m.id == id; };
        const bool alive = std::any_of(plots_.begin(), plots_.end(), same) || std::any_of(dets_.begin(), dets_.end(), same);
        it = alive ? std::next(it) : pinned_.erase(it);
    }
    update();
}

/* Default: a plot disappears when the antenna has turned a full revolution
 * since it pointed at the plot, i.e. it sweeps over the plot azimuth in the
 * next scan (safety net for gaps in the data: two scans later it is gone in
 * any case). Options: two revolutions, or all plots until cleared. */
bool PpiWidget::expired(const PlotMark &m) const
{
    if (history_ == History::All) return false;
    const double turns = history_ == History::TwoScans ? 2.0 : 1.0;
    const long long scans = (long long)antenna_.scan - (long long)m.scan;
    return turned(m) >= 360.0 * turns || scans >= 1 + (long long)turns || scans <= -2;
}

bool PpiWidget::shown(const PlotMark &m) const
{
    const bool modeOk = m.isDet ? mode_ != DisplayMode::Plots : mode_ != DisplayMode::Detections;
    return modeOk && (m.pulse == PULSE_LONG ? showLong_ : showShort_) && !(m.det.snr_dB < minSnr_dB_);
}

size_t PpiWidget::visibleCount() const
{
    return (size_t)std::count_if(plots_.begin(), plots_.end(), [this](const PlotMark &m) { return shown(m); });
}

size_t PpiWidget::visibleDetectionCount() const
{
    return (size_t)std::count_if(dets_.begin(), dets_.end(), [this](const PlotMark &m) { return shown(m); });
}

const PlotMark *PpiWidget::pick(QPointF screen) const
{
    const PlotMark *best = nullptr;
    double bestD = kPickPx;
    for (const std::deque<PlotMark> *list : { &plots_, &dets_ }) {   /* plots first */
        for (const PlotMark &m : *list) {
            if (!shown(m)) continue;
            const QPointF d = toScreen(m.east_m, m.north_m) - screen;
            const double dist = std::hypot(d.x(), d.y());
            if (dist < bestD || (dist <= bestD && !best)) { bestD = dist; best = &m; }
        }
    }
    return best;
}

void PpiWidget::clickAt(QPointF screen)
{
    const PlotMark *m = pick(screen);
    if (!m) return;
    const uint64_t id = m->id;
    if (!pinned_.erase(id)) pinned_.insert(id);
    update();
}

QString PpiWidget::saveSnapshotDialog()
{
    const QString name = QStringLiteral("ppi_%1.png").arg(QDateTime::currentDateTime().toString(QStringLiteral("yyyyMMdd_HHmmss")));
    const QString path = QFileDialog::getSaveFileName(this, QStringLiteral("Save a screenshot of the PPI"), name,
                                                      QStringLiteral("PNG (*.png)"));
    if (path.isEmpty() || !saveSnapshot(path)) return QString();
    return path;
}

bool PpiWidget::saveSnapshot(const QString &path)
{
    const QPixmap pixmap = grab();
    if (pixmap.save(path, "PNG")) return true;
    std::fprintf(stderr, "PPI: cannot save the snapshot to %s\n", path.toLocal8Bit().constData());
    return false;
}

/* ---- view ---------------------------------------------------------------- */

void PpiWidget::setView(double zoom, double e, double n)
{
    zoom_ = std::max(1.0, std::min(kMaxZoom, zoom));
    /* keep the instrumented range in view: the view centre stays within it */
    const double r = std::hypot(e, n), rMax = maxRange_m_ * (1.0 - 1.0 / zoom_);
    if (r > rMax) {
        e *= rMax / r;
        n *= rMax / r;
    }
    viewE_ = e;
    viewN_ = n;
    backgroundDirty_ = true;
    update();
}

void PpiWidget::zoomAt(double factor, QPointF screen)
{
    layoutDisplay();
    const QPointF g = toGround(screen);          /* this ground point stays under the cursor */
    const double z = std::max(1.0, std::min(kMaxZoom, zoom_ * factor));
    const double k = radius_px_ * z / maxRange_m_;
    setView(z, g.x() - (screen.x() - centre_.x()) / k, g.y() + (screen.y() - centre_.y()) / k);
}

void PpiWidget::panPixels(QPointF delta)
{
    layoutDisplay();
    const double k = scale();
    setView(zoom_, viewE_ - delta.x() / k, viewN_ + delta.y() / k);
}

void PpiWidget::resetView()
{
    setView(1.0, 0.0, 0.0);
}

/* ---- events -------------------------------------------------------------- */

void PpiWidget::closeEvent(QCloseEvent *event)
{
    if (!finished_ && stopFlag_) *stopFlag_ = 1;   /* ask the worker to stop */
    event->accept();
}

void PpiWidget::wheelEvent(QWheelEvent *event)
{
    const double steps = event->angleDelta().y() / 120.0;
    if (steps != 0.0) zoomAt(std::pow(1.25, steps), event_pos(event));
    event->accept();
}

void PpiWidget::mousePressEvent(QMouseEvent *event)
{
    if (event->button() != Qt::LeftButton) return;
    dragging_ = true;
    dragMoved_ = false;
    dragStart_ = event_pos(event);
    dragE_ = viewE_;
    dragN_ = viewN_;
}

void PpiWidget::mouseMoveEvent(QMouseEvent *event)
{
    const QPointF pos = event_pos(event);
    if (dragging_) {
        const QPointF d = pos - dragStart_;
        if (std::fabs(d.x()) > 3.0 || std::fabs(d.y()) > 3.0) dragMoved_ = true;
        if (dragMoved_) {
            const double k = scale();
            setView(zoom_, dragE_ - d.x() / k, dragN_ + d.y() / k);
            setCursor(Qt::ClosedHandCursor);
        }
    }
    updateCursor(pos);
    uint64_t id = 0;
    const QString text = dragging_ ? QString() : hoverText(pos, &id);
    if (id != hoverId_) {
        hoverId_ = id;
        if (id) QToolTip::showText(mapToGlobal(pos.toPoint()), text, this);
        else QToolTip::hideText();
        update();
    }
}

void PpiWidget::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() != Qt::LeftButton || !dragging_) return;
    if (!dragMoved_) clickAt(event_pos(event));
    dragging_ = false;
    unsetCursor();
}

void PpiWidget::mouseDoubleClickEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) resetView();
}

void PpiWidget::leaveEvent(QEvent *event)
{
    (void)event;
    cursorText_.clear();
    hoverId_ = 0;
    update();
}

void PpiWidget::keyPressEvent(QKeyEvent *event)
{
    const QPointF mid = centre_;
    const double step = 0.1 * std::min(width(), height());
    switch (event->key()) {
    case Qt::Key_Plus: case Qt::Key_Equal: zoomAt(1.5, mid); break;
    case Qt::Key_Minus:                    zoomAt(1.0 / 1.5, mid); break;
    case Qt::Key_Left:  panPixels(QPointF(step, 0.0)); break;
    case Qt::Key_Right: panPixels(QPointF(-step, 0.0)); break;
    case Qt::Key_Up:    panPixels(QPointF(0.0, step)); break;
    case Qt::Key_Down:  panPixels(QPointF(0.0, -step)); break;
    case Qt::Key_R: case Qt::Key_Home: resetView(); break;
    case Qt::Key_C:     pinned_.clear(); update(); break;
    case Qt::Key_T:     setShowTracks(!showTracks_); break;
    case Qt::Key_D:     setDisplayMode(mode_ == DisplayMode::Plots ? DisplayMode::Detections
                                       : mode_ == DisplayMode::Detections ? DisplayMode::Both : DisplayMode::Plots);
                        break;
    case Qt::Key_Space: setPaused(!paused_); break;
    case Qt::Key_S: {
        saveSnapshotDialog();
        break;
    }
    default: QWidget::keyPressEvent(event);
    }
}

void PpiWidget::contextMenuEvent(QContextMenuEvent *event)
{
    QMenu menu(this);
    auto toggle = [&menu](const QString &text, bool on, std::function<void(bool)> f) {
        QAction *a = menu.addAction(text);
        a->setCheckable(true);
        a->setChecked(on);
        QObject::connect(a, &QAction::toggled, f);
    };
    QMenu *modeMenu = menu.addMenu(QStringLiteral("Show (D)"));
    QActionGroup *modeGroup = new QActionGroup(modeMenu);
    const std::pair<const char *, DisplayMode> modeItems[] = {
        { "Plots (after the plot extractor)", DisplayMode::Plots },
        { "CFAR detections (before the extractor)", DisplayMode::Detections },
        { "Plots and detections", DisplayMode::Both } };
    for (const auto &item : modeItems) {
        QAction *a = modeMenu->addAction(QString::fromLatin1(item.first));
        a->setCheckable(true);
        a->setChecked(mode_ == item.second);
        modeGroup->addAction(a);
        const DisplayMode m = item.second;
        QObject::connect(a, &QAction::triggered, [this, m]() { setDisplayMode(m); });
    }
    toggle(QStringLiteral("Short pulse"), showShort_, [this](bool on) { setShowShort(on); });
    toggle(QStringLiteral("Long pulse"), showLong_, [this](bool on) { setShowLong(on); });
    toggle(QStringLiteral("Plot labels (range, velocity)"), showLabels_, [this](bool on) { showLabels_ = on; update(); });
    toggle(QStringLiteral("Range rings and spokes"), showRings_,
           [this](bool on) { showRings_ = on; backgroundDirty_ = true; update(); });
    toggle(QStringLiteral("Antenna trail"), showTrail_, [this](bool on) { showTrail_ = on; update(); });
    toggle(QStringLiteral("Tracks (T)"), showTracks_, [this](bool on) { setShowTracks(on); });

    QMenu *hist = menu.addMenu(QStringLiteral("Plots shown"));
    QActionGroup *histGroup = new QActionGroup(hist);
    const std::pair<const char *, History> histItems[] = {
        { "Until the antenna passes again", History::UntilSwept },
        { "Two scans", History::TwoScans },
        { "All (until cleared)", History::All } };
    for (const auto &item : histItems) {
        QAction *a = hist->addAction(QString::fromLatin1(item.first));
        a->setCheckable(true);
        a->setChecked(history_ == item.second);
        histGroup->addAction(a);
        const History h = item.second;
        QObject::connect(a, &QAction::triggered, [this, h]() { setHistory(h); });
    }
    QMenu *snr = menu.addMenu(QStringLiteral("Minimum SNR"));
    QActionGroup *snrGroup = new QActionGroup(snr);
    for (double dB : { -1e9, 10.0, 13.0, 15.0, 20.0, 25.0, 30.0 }) {
        QAction *a = snr->addAction(dB < -1e8 ? QStringLiteral("All") : QStringLiteral("%1 dB").arg(dB));
        a->setCheckable(true);
        a->setChecked(minSnr_dB_ == dB);
        snrGroup->addAction(a);
        QObject::connect(a, &QAction::triggered, [this, dB]() { setMinSnr(dB); });
    }
    menu.addSeparator();
    QObject::connect(menu.addAction(QStringLiteral("Reset view (double click)")), &QAction::triggered,
                     [this]() { resetView(); });
    QObject::connect(menu.addAction(QStringLiteral("Clear display")), &QAction::triggered, [this]() { clearPlots(); });
    QObject::connect(menu.addAction(QStringLiteral("Remove data tips (C)")), &QAction::triggered,
                     [this]() { pinned_.clear(); update(); });
    QObject::connect(menu.addAction(paused_ ? QStringLiteral("Resume display (Space)") : QStringLiteral("Pause display (Space)")),
                     &QAction::triggered, [this]() { setPaused(!paused_); });
    QObject::connect(menu.addAction(QStringLiteral("Save screenshot ... (S)")), &QAction::triggered,
                     [this]() { saveSnapshotDialog(); });
    menu.exec(event->globalPos());
}

void PpiWidget::updateCursor(QPointF screen)
{
    if (!viewRect_.contains(screen)) {
        cursorText_.clear();
    } else {
        const QPointF g = toGround(screen);
        cursorText_ = QStringLiteral("Cursor %1, %2   (E %3, N %4)")
                          .arg(km_text(std::hypot(g.x(), g.y())))
                          .arg(deg_text(wrap360(std::atan2(g.x(), g.y()) / kDeg), 1))
                          .arg(km_text(g.x()))
                          .arg(km_text(g.y()));
    }
    update();
}

/* ---- drawing ------------------------------------------------------------- */

/* Square polar display below the status text, with room for the labels. */
void PpiWidget::layoutDisplay()
{
    const double w = width(), h = height();
    const double top = kHeaderPx + kLabelPx;
    const double side = std::max(50.0, std::min(w - 2.0 * kLabelPx, h - top - kLabelPx));
    radius_px_ = side / 2.0;
    centre_ = QPointF(w / 2.0, top + (h - top - kLabelPx) / 2.0);
    viewRect_ = QRectF(0.0, kHeaderPx, w, h - kHeaderPx);
}

/* Farthest point of the ray from the radar along az_deg (up to rEnd_px) that
 * lies inside rect; false if the ray misses the rectangle. */
bool ray_end_in_rect(QPointF origin, double az_deg, double rEnd_px, const QRectF &rect, QPointF *end)
{
    const double dx = std::sin(az_deg * kDeg), dy = -std::cos(az_deg * kDeg);
    double t0 = 0.0, t1 = rEnd_px;                         /* Liang-Barsky clipping */
    const double p[4] = { -dx, dx, -dy, dy };
    const double q[4] = { origin.x() - rect.left(), rect.right() - origin.x(),
                          origin.y() - rect.top(), rect.bottom() - origin.y() };
    for (int i = 0; i < 4; ++i) {
        if (std::fabs(p[i]) < 1e-12) {
            if (q[i] < 0.0) return false;
            continue;
        }
        const double t = q[i] / p[i];
        if (p[i] < 0.0) t0 = std::max(t0, t);
        else t1 = std::min(t1, t);
    }
    if (t0 > t1) return false;
    *end = QPointF(origin.x() + t1 * dx, origin.y() + t1 * dy);
    return true;
}

void PpiWidget::drawBackground()
{
    const qreal dpr = devicePixelRatioF();
    background_ = QPixmap(size() * dpr);
    background_.setDevicePixelRatio(dpr);
    background_.fill(Qt::black);
    backgroundSize_ = size();
    backgroundDirty_ = false;

    QPainter p(&background_);
    p.setRenderHint(QPainter::Antialiasing, true);
    p.setClipRect(viewRect_);
    const QFont font = pixel_font(this->font(), 12);
    const QFontMetrics fm(font);
    p.setFont(font);
    const double k = scale(), rim_px = maxRange_m_ * k;
    const QPointF origin = toScreen(0.0, 0.0);

    if (showRings_) {
        /* azimuth spokes every 30 deg; labels outside the rim or, when zoomed,
         * at the edge of the view */
        const QRectF inner = viewRect_.adjusted(14.0, 10.0, -14.0, -10.0);
        for (int az = 0; az < 360; az += 30) {
            QPointF end;
            if (!ray_end_in_rect(origin, az, rim_px, viewRect_, &end)) continue;
            p.setPen(QPen(kColorGrid, 1.0));
            p.drawLine(origin, polar(az, maxRange_m_));
            QPointF c;
            if (zoom_ <= 1.0 + 1e-9) c = polar(az, maxRange_m_ * (1.0 + 15.0 / rim_px));
            else if (!ray_end_in_rect(origin, az, rim_px, inner, &c)) continue;
            const QString text = deg_text(az, 0);
            p.setPen(kColorLabel);
            p.drawText(QPointF(c.x() - fm.horizontalAdvance(text) / 2.0, c.y() + fm.ascent() / 2.0 - 1.0), text);
        }
        /* range rings for the visible range, labels between the 0 and 30 deg spokes */
        const double step = ring_step_m(maxRange_m_ / zoom_);
        for (double r = step; r <= maxRange_m_ * (1.0 + 1e-9); r += step) {
            p.setPen(QPen(kColorGrid, 1.0));
            p.setBrush(Qt::NoBrush);
            p.drawEllipse(origin, r * k, r * k);
            p.setPen(kColorLabel);
            p.drawText(polar(15.0, r) + QPointF(3.0, -3.0), km_text(r));
        }
    }
    /* rim of the instrumented range */
    p.setPen(QPen(kColorRim, 1.5));
    p.setBrush(Qt::NoBrush);
    p.drawEllipse(origin, rim_px, rim_px);
}

/* Antenna line with a fading trail behind it (the antenna turns clockwise). */
void PpiWidget::drawAntenna(QPainter &p) const
{
    if (!antenna_.valid) return;
    const QPointF origin = toScreen(0.0, 0.0);
    const double rim_px = maxRange_m_ * scale();
    if (showTrail_) {
        const QRectF disc(origin.x() - rim_px, origin.y() - rim_px, 2.0 * rim_px, 2.0 * rim_px);
        const double qtAngle = 90.0 - antenna_.az_deg;   /* Qt: counter-clockwise from East */
        QConicalGradient trail(origin, qtAngle);
        QColor head(kColorAntenna), tail(kColorAntenna);
        head.setAlpha(90);
        tail.setAlpha(0);
        trail.setColorAt(0.0, head);
        trail.setColorAt(kTrailDeg / 360.0, tail);
        trail.setColorAt(1.0, tail);
        p.setPen(Qt::NoPen);
        p.setBrush(trail);
        p.drawPie(disc, (int)std::lround(qtAngle * 16.0), (int)std::lround(kTrailDeg * 16.0));
    }
    p.setPen(QPen(kColorAntenna, 2.0));
    p.drawLine(origin, polar(antenna_.az_deg, maxRange_m_));
}

/* Filled circles; older plots are dimmer (phosphor-like decay over a scan). */
void PpiWidget::drawPlots(QPainter &p) const
{
    const QFont font = pixel_font(this->font(), 11);
    p.setFont(font);
    for (const PlotMark &m : plots_) {
        if (!shown(m)) continue;
        const QPointF s = toScreen(m.east_m, m.north_m);
        if (!viewRect_.adjusted(-10, -10, 10, 10).contains(s)) continue;
        QColor c = (m.pulse == PULSE_LONG) ? kColorLong : kColorShort;
        const double age = std::min(1.0, turned(m) / 360.0);
        c.setAlpha(255 - (int)(130.0 * age));
        p.setPen(Qt::NoPen);
        p.setBrush(c);
        p.drawEllipse(s, kPlotRadius, kPlotRadius);
        if (m.id == hoverId_) {                         /* highlight the plot under the cursor */
            p.setPen(QPen(Qt::white, 1.5));
            p.setBrush(Qt::NoBrush);
            p.drawEllipse(s, kPlotRadius + 4.0, kPlotRadius + 4.0);
        }
        if (showLabels_) {
            p.setPen(kColorLabel);
            p.drawText(s + QPointF(7.0, 4.0), QStringLiteral("%1 km, %2 m/s")
                                                  .arg(m.det.range_m / 1e3, 0, 'f', 1)
                                                  .arg(m.det.vr_mps, 0, 'f', 1));
        }
    }
}

/* Hollow triangle around the position of the last update (the plot lies
 * inside it), velocity leader and number; tentative tracks in gray. */
void PpiWidget::drawTracks(QPainter &p) const
{
    if (!showTracks_) return;
    p.setFont(pixel_font(this->font(), 11));
    for (const auto &e : tracks_) {
        const Track &t = e.second;
        const QPointF s = toScreen(t.pos_m[0], t.pos_m[1]);
        if (!viewRect_.adjusted(-40, -40, 40, 40).contains(s)) continue;
        const QColor c = t.confirmed ? kColorTrack : kColorTentative;
        const double h = kTrackPx;
        const QPointF tri[3] = { s + QPointF(0.0, -h), s + QPointF(h * 0.866, h * 0.5), s + QPointF(-h * 0.866, h * 0.5) };
        p.setPen(QPen(c, t.confirmed ? 1.6 : 1.0));
        p.setBrush(Qt::NoBrush);
        p.drawPolygon(tri, 3);
        if (!t.confirmed) continue;
        p.drawLine(s, toScreen(t.pos_m[0] + kLeader_s * t.vel_mps[0], t.pos_m[1] + kLeader_s * t.vel_mps[1]));
        p.drawText(s + QPointF(h, h + 4.0), QStringLiteral("T%1").arg(t.id));
        if (hoverId_ == ((1ull << 63) | t.id)) {
            p.setPen(QPen(Qt::white, 1.5));
            p.drawEllipse(s, h + 4.0, h + 4.0);
        }
    }
}

size_t PpiWidget::confirmedCount() const
{
    return (size_t)std::count_if(tracks_.begin(), tracks_.end(),
                                 [](const std::pair<const uint32_t, Track> &e) { return e.second.confirmed != 0; });
}

const Track *PpiWidget::pickTrack(QPointF screen, double *dist) const
{
    if (!showTracks_) return nullptr;
    const Track *best = nullptr;
    double bestD = kPickPx;
    for (const auto &e : tracks_) {
        const QPointF d = toScreen(e.second.pos_m[0], e.second.pos_m[1]) - screen;
        const double r = std::hypot(d.x(), d.y());
        if (r <= bestD) { bestD = r; best = &e.second; }
    }
    if (dist) *dist = bestD;
    return best;
}

/* Data tip of a track. */
QString track_text(const Track &t)
{
    QStringList rows;
    const double R = std::sqrt(t.pos_m[0] * t.pos_m[0] + t.pos_m[1] * t.pos_m[1] + t.pos_m[2] * t.pos_m[2]);
    rows << QStringLiteral("<b>Track %1 (%2, %3-D EKF)</b>").arg(t.id).arg(t.confirmed ? "confirmed" : "tentative").arg(t.dim)
         << QStringLiteral("Range: %1 km").arg(R / 1e3, 0, 'f', 3)
         << QStringLiteral("Azimuth: %1").arg(deg_text(wrap360(std::atan2(t.pos_m[0], t.pos_m[1]) / kDeg), 2));
    if (t.dim == 3)
        rows << QStringLiteral("Height: %1 m").arg(t.alt_m, 0, 'f', 0)
             << QStringLiteral("Climb: %1 m/s").arg(t.vel_mps[2], 0, 'f', 1);
    rows << QStringLiteral("Speed: %1 m/s").arg(std::hypot(t.vel_mps[0], t.vel_mps[1]), 0, 'f', 1)
         << QStringLiteral("Heading: %1").arg(deg_text(wrap360(std::atan2(t.vel_mps[0], t.vel_mps[1]) / kDeg), 1));
    const double acc = std::hypot(t.acc_mps2[0], t.acc_mps2[1]);
    if (acc > 0.0) rows << QStringLiteral("Acceleration: %1 m/s\u00b2").arg(acc, 0, 'f', 2);   /* CA model */
    rows << QStringLiteral("Updates: %1 (missed scans %2)").arg(t.hits).arg(t.misses)
         << QStringLiteral("Last update: %1 s").arg(t.time_s, 0, 'f', 3);
    return rows.join(QStringLiteral("<br>"));
}

/* Data tip under the cursor and a key for it (plot/detection id, or the track
 * number with the top bit set; 0 = nothing). */
QString PpiWidget::hoverText(QPointF screen, uint64_t *key) const
{
    double dTrack = 0.0;
    const Track *t = pickTrack(screen, &dTrack);
    const PlotMark *m = pick(screen);
    if (m && t) {
        const QPointF d = toScreen(m->east_m, m->north_m) - screen;
        if (std::hypot(d.x(), d.y()) > dTrack - 3.0) m = nullptr;   /* the track wins */
    }
    if (m) { *key = m->id; return plot_text(m->det, m->isDet); }
    if (t) { *key = (1ull << 63) | t->id; return track_text(*t); }
    *key = 0;
    return QString();
}

/* CFAR detections: an arc over the azimuth step of the CPI update at the
 * ground range of the detection (a dot when the arc is shorter than a pixel). */
void PpiWidget::drawDetections(QPainter &p) const
{
    if (mode_ == DisplayMode::Plots) return;
    const QPointF origin = toScreen(0.0, 0.0);
    const double k = scale();
    const QRectF visible = viewRect_.adjusted(-10, -10, 10, 10);
    for (const PlotMark &m : dets_) {
        if (!shown(m)) continue;
        const QPointF s = toScreen(m.east_m, m.north_m);
        const double g_px = std::hypot(m.east_m, m.north_m) * k, arc_px = g_px * m.arc_deg * kDeg;
        if (!visible.adjusted(-arc_px, -arc_px, arc_px, arc_px).contains(s)) continue;
        QColor c = (m.pulse == PULSE_LONG) ? kColorDetLong : kColorDetShort;
        c.setAlpha(255 - (int)(130.0 * std::min(1.0, turned(m) / 360.0)));
        p.setPen(QPen(c, 3.0, Qt::SolidLine, Qt::RoundCap));
        p.setBrush(Qt::NoBrush);
        if (arc_px < 1.5) {
            p.drawPoint(s);
        } else {                                          /* Qt angles: counter-clockwise from East */
            const QRectF circle(origin.x() - g_px, origin.y() - g_px, 2.0 * g_px, 2.0 * g_px);
            const double az = std::atan2(m.east_m, m.north_m) / kDeg;
            QPainterPath path;
            path.arcMoveTo(circle, 90.0 - az - m.arc_deg / 2.0);
            path.arcTo(circle, 90.0 - az - m.arc_deg / 2.0, m.arc_deg);
            p.drawPath(path);
        }
        if (m.id == hoverId_) {
            p.setPen(QPen(Qt::white, 1.5));
            p.drawEllipse(s, kPlotRadius + 4.0, kPlotRadius + 4.0);
        }
    }
}

/* Pinned data tips: a box with the values next to the plot or detection. */
void PpiWidget::drawPins(QPainter &p) const
{
    const QFont font = pixel_font(this->font(), 11);
    const QFontMetrics fm(font);
    p.setFont(font);
    std::vector<const PlotMark *> marks;
    for (const std::deque<PlotMark> *list : { &plots_, &dets_ })
        for (const PlotMark &mk : *list)
            if (pinned_.count(mk.id)) marks.push_back(&mk);
    for (const PlotMark *pm : marks) {
        const PlotMark &m = *pm;
        const QPointF s = toScreen(m.east_m, m.north_m);
        QString text = plot_text(m.det, m.isDet);
        text.replace(QStringLiteral("<b>"), QString()).replace(QStringLiteral("</b>"), QString());
        const QStringList rows = text.split(QStringLiteral("<br>"));
        int w = 0;
        for (const QString &r : rows) w = std::max(w, fm.horizontalAdvance(r));
        const double rowH = fm.height(), boxW = w + 12.0, boxH = rows.size() * rowH + 8.0;
        QRectF box(s.x() + 12.0, s.y() + 8.0, boxW, boxH);
        if (box.right() > width()) box.moveRight(s.x() - 12.0);          /* stay inside the widget */
        if (box.bottom() > height()) box.moveBottom(s.y() - 8.0);
        p.setPen(QPen(kColorPin, 1.5));
        p.setBrush(Qt::NoBrush);
        p.drawRect(QRectF(s.x() - 7.0, s.y() - 7.0, 14.0, 14.0));
        p.drawLine(s, box.topLeft());
        p.setBrush(QColor(10, 40, 10, 225));
        p.drawRect(box);
        p.setPen(Qt::white);
        for (int i = 0; i < rows.size(); ++i)
            p.drawText(QPointF(box.left() + 6.0, box.top() + 4.0 + fm.ascent() + i * rowH), rows[i]);
    }
}

void PpiWidget::drawStatus(QPainter &p) const
{
    QString line1;
    if (antenna_.valid) {
        line1 = QStringLiteral("Scan %1    t = %2 s    Antenna az %3")
                    .arg(antenna_.scan)
                    .arg(antenna_.time_s, 0, 'f', 2)
                    .arg(deg_text(antenna_.az_deg, 1));
        if (mode_ != DisplayMode::Detections)
            line1 += QStringLiteral("    Plots %1 / %2").arg((qulonglong)visibleCount()).arg((qulonglong)plots_.size());
        if (mode_ != DisplayMode::Plots)
            line1 += QStringLiteral("    Detections %1").arg((qulonglong)visibleDetectionCount());
        if (showTracks_ && !tracks_.empty())
            line1 += QStringLiteral("    Tracks %1").arg((qulonglong)confirmedCount());
    }
    else
        line1 = finished_ ? QStringLiteral("No antenna data received")
                          : QStringLiteral("Waiting for data ...");
    if (paused_) line1 += QStringLiteral("    [PAUSED]");

    QString line2 = QStringLiteral("Range %1, zoom x%2, rings every %3    ")
                        .arg(km_text(maxRange_m_))
                        .arg(zoom_, 0, 'f', 1)
                        .arg(km_text(ring_step_m(maxRange_m_ / zoom_)));
    if (!cursorText_.isEmpty())
        line2 += cursorText_;
    else if (!finished_)
        line2 += QStringLiteral("receiving");
    else if (workerCode_ == 0)
        line2 += QStringLiteral("processing finished");
    else
        line2 += QStringLiteral("processing stopped with error %1").arg(workerCode_);

    p.setPen(Qt::white);
    p.setFont(pixel_font(font(), 16, true));
    p.drawText(QPointF(12.0, 22.0), line1);
    p.setPen(kColorLabel);
    p.setFont(pixel_font(font(), 12));
    p.drawText(QPointF(12.0, 41.0), line2);
}

/* Small legend in the bottom-left corner, with the mouse help. */
void PpiWidget::drawLegend(QPainter &p) const
{
    enum Kind { PlotS, PlotL, DetS, DetL, TrackK, Antenna, Help };
    std::vector<std::pair<Kind, QString>> rows;
    if (mode_ != DisplayMode::Detections) {
        rows.emplace_back(PlotS, QStringLiteral("Plot, short pulse"));
        rows.emplace_back(PlotL, QStringLiteral("Plot, long pulse"));
    }
    if (mode_ != DisplayMode::Plots) {
        rows.emplace_back(DetS, QStringLiteral("CFAR detection, short pulse"));
        rows.emplace_back(DetL, QStringLiteral("CFAR detection, long pulse"));
    }
    if (showTracks_) rows.emplace_back(TrackK, QStringLiteral("Track (gray: tentative)"));
    rows.emplace_back(Antenna, QStringLiteral("Antenna"));
    rows.emplace_back(Help, QStringLiteral("Wheel zoom, drag pan, right click options"));

    const QFont font = pixel_font(this->font(), 12);
    const QFontMetrics fm(font);
    int textWidth = 0;
    for (const auto &r : rows) textWidth = std::max(textWidth, fm.horizontalAdvance(r.second));
    const double rowH = 18.0, pad = 8.0, n = (double)rows.size();
    const QRectF box(8.0, height() - 8.0 - n * rowH - 2.0 * pad, 2.0 * pad + 26.0 + textWidth, n * rowH + 2.0 * pad);
    p.setPen(QPen(kColorGrid, 1.0));
    p.setBrush(QColor(20, 20, 20, 220));
    p.drawRect(box);

    p.setFont(font);
    for (size_t i = 0; i < rows.size(); ++i) {
        const Kind kind = rows[i].first;
        const double y = box.top() + pad + rowH * ((double)i + 0.5);
        const QPointF mark(box.left() + pad + 8.0, y);
        if (kind == PlotS || kind == PlotL) {
            p.setPen(Qt::NoPen);
            p.setBrush(kind == PlotS ? kColorShort : kColorLong);
            p.drawEllipse(mark, kPlotRadius, kPlotRadius);
        } else if (kind == DetS || kind == DetL) {
            p.setPen(QPen(kind == DetS ? kColorDetShort : kColorDetLong, 3.0, Qt::SolidLine, Qt::RoundCap));
            p.setBrush(Qt::NoBrush);
            p.drawArc(QRectF(mark.x() - 8.0, mark.y() - 2.0, 16.0, 16.0), 50 * 16, 80 * 16);
        } else if (kind == TrackK) {
            const double h = 7.0;
            const QPointF tri[3] = { mark + QPointF(0.0, -h), mark + QPointF(h * 0.866, h * 0.5),
                                     mark + QPointF(-h * 0.866, h * 0.5) };
            p.setPen(QPen(kColorTrack, 1.5));
            p.setBrush(Qt::NoBrush);
            p.drawPolygon(tri, 3);
        } else if (kind == Antenna) {
            p.setPen(QPen(kColorAntenna, 2.0));
            p.drawLine(mark - QPointF(8.0, 0.0), mark + QPointF(8.0, 0.0));
        }
        p.setPen(kind == Help ? kColorLabel : QColor(Qt::white));
        p.drawText(QPointF(box.left() + pad + (kind == Help ? 0.0 : 26.0), y + fm.ascent() / 2.0 - 1.0), rows[i].second);
    }
}

void PpiWidget::paintEvent(QPaintEvent *event)
{
    (void)event;
    const double t0 = gui_seconds();
    struct AddTime {                    /* also when a draw call returns early */
        PpiWidget *w;
        double t0;
        ~AddTime() { w->addGuiMs((gui_seconds() - t0) * 1e3); }
    } addTime{ this, t0 };
    layoutDisplay();
    if (background_.isNull() || backgroundDirty_ || backgroundSize_ != size()) drawBackground();

    QPainter p(this);
    p.drawPixmap(0, 0, background_);
    p.setRenderHint(QPainter::Antialiasing, true);
    p.save();
    p.setClipRect(viewRect_);
    drawAntenna(p);
    drawDetections(p);
    drawPlots(p);
    drawTracks(p);
    p.restore();
    drawPins(p);
    drawStatus(p);
    drawLegend(p);
}

bool save_snapshot(PpiWidget &widget, const char *path)
{
    return widget.saveSnapshot(QString::fromLocal8Bit(path));
}

/* ------------------------------------------------------------------------ */
/* Window: the PPI and a panel with the display mode, the screenshot button  */
/* and the CFAR threshold                                                    */
/* ------------------------------------------------------------------------ */

class PpiWindow : public QWidget {
public:
    PpiWindow(double maxRange_m, volatile sig_atomic_t *stopFlag) : stopFlag_(stopFlag)
    {
        ppi_ = new PpiWidget(maxRange_m, nullptr);       /* the window handles closing */
        QWidget *panel = new QWidget(this);
        panel->setFixedWidth(230);
        panel->setStyleSheet(QStringLiteral(
            "QWidget { background: #181818; color: #e0e0e0; }"
            "QPushButton, QComboBox, QDoubleSpinBox { background: #2a2a2a; border: 1px solid #555; padding: 4px; }"
            "QPushButton:hover { background: #383838; }"
            "QLabel#title { font-weight: bold; color: white; }"
            "QLabel#info { color: #a0e0a0; }"));
        QVBoxLayout *v = new QVBoxLayout(panel);
        auto title = [panel, v](const QString &text) {
            QLabel *l = new QLabel(text, panel);
            l->setObjectName(QStringLiteral("title"));
            v->addSpacing(8);
            v->addWidget(l);
        };

        title(QStringLiteral("Display"));
        mode_ = new QComboBox(panel);
        mode_->addItem(QStringLiteral("Plots (after extractor)"));
        mode_->addItem(QStringLiteral("Detections (before extractor)"));
        mode_->addItem(QStringLiteral("Plots and detections"));
        v->addWidget(mode_);
        QObject::connect(mode_, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
                         [this](int i) { ppi_->setDisplayMode((DisplayMode)i); ppi_->setFocus(); });
        tracks_ = new QCheckBox(QStringLiteral("Tracks (T)"), panel);
        tracks_->setChecked(ppi_->showTracks());
        v->addWidget(tracks_);
        QObject::connect(tracks_, &QCheckBox::toggled, [this](bool on) { ppi_->setShowTracks(on); ppi_->setFocus(); });

        title(QStringLiteral("Screenshot"));
        save_ = new QPushButton(QStringLiteral("Save screenshot (PNG) ..."), panel);
        v->addWidget(save_);
        saved_ = new QLabel(panel);
        saved_->setObjectName(QStringLiteral("info"));
        saved_->setWordWrap(true);
        v->addWidget(saved_);
        QObject::connect(save_, &QPushButton::clicked, [this]() {
            const QString path = ppi_->saveSnapshotDialog();
            if (!path.isEmpty()) saved_->setText(QStringLiteral("saved: %1").arg(path));
        });

        title(QStringLiteral("CFAR threshold [dB]"));
        cfar_ = new QDoubleSpinBox(panel);
        cfar_->setRange(0.0, 50.0);
        cfar_->setDecimals(2);
        cfar_->setSingleStep(0.5);
        cfar_->setEnabled(false);
        v->addWidget(cfar_);
        QHBoxLayout *h = new QHBoxLayout();
        apply_ = new QPushButton(QStringLiteral("Apply"), panel);
        default_ = new QPushButton(QStringLiteral("Default"), panel);
        apply_->setEnabled(false);
        default_->setEnabled(false);
        h->addWidget(apply_);
        h->addWidget(default_);
        v->addLayout(h);
        cfarInfo_ = new QLabel(QStringLiteral("waiting for the processing ..."), panel);
        cfarInfo_->setObjectName(QStringLiteral("info"));
        cfarInfo_->setWordWrap(true);
        cfarInfo_->setMinimumHeight(cfarInfo_->fontMetrics().lineSpacing() * 2 + 4);
        v->addWidget(cfarInfo_);
        QObject::connect(apply_, &QPushButton::clicked, [this]() { requestCfar(cfar_->value()); });
        QObject::connect(default_, &QPushButton::clicked, [this]() { requestCfar(cfarDefault_dB_); });

        v->addStretch(1);
        title(QStringLiteral("Mean time per PRI [us]"));
        timing_ = new QLabel(QStringLiteral("waiting for the processing ..."), panel);
        QFont mono = QFontDatabase::systemFont(QFontDatabase::FixedFont);
        mono.setPixelSize(11);
        timing_->setFont(mono);
        timing_->setStyleSheet(QStringLiteral("background: #101010; padding: 4px;"));
        timing_->setMinimumHeight(QFontMetrics(mono).lineSpacing() * 12 + 12);  /* 12 rows, stable layout */
        timing_->setAlignment(Qt::AlignLeft | Qt::AlignTop);
        v->addWidget(timing_);
        QLabel *help = new QLabel(QStringLiteral("Wheel: zoom   Drag: pan\nDouble click: reset\nHover: data tip   "
                                                 "Click: pin\nKeys: + - arrows R C D T S Space\nRight click: more options"),
                                  panel);
        help->setStyleSheet(QStringLiteral("color: #909090;"));
        v->addWidget(help);

        QHBoxLayout *top = new QHBoxLayout(this);
        top->setContentsMargins(0, 0, 0, 0);
        top->setSpacing(0);
        top->addWidget(ppi_, 1);
        top->addWidget(panel);
        ppi_->setFocus();
    }

    PpiWidget &ppi() { return *ppi_; }
    QCheckBox *tracksBox() { return tracks_; }
    QString timingText() const { return timing_->text(); }
    QComboBox *modeBox() { return mode_; }
    QPushButton *applyButton() { return apply_; }
    QPushButton *defaultButton() { return default_; }
    QDoubleSpinBox *cfarBox() { return cfar_; }
    QString cfarInfo() const { return cfarInfo_->text(); }

    void setFinished(int workerCode)
    {
        finished_ = true;
        ppi_->setFinished(workerCode);
    }

    /* GUI timer: new data for the PPI, threshold reported by the worker */
    void refresh()
    {
        const double t0 = gui_seconds();
        ppi_->refresh();
        ppi_->addGuiMs((gui_seconds() - t0) * 1e3);
        if (tracks_->isChecked() != ppi_->showTracks()) {          /* changed with T or the context menu */
            const bool b = tracks_->blockSignals(true);
            tracks_->setChecked(ppi_->showTracks());
            tracks_->blockSignals(b);
        }
        showTiming();
        if (mode_->currentIndex() != (int)ppi_->displayMode()) {   /* changed with D or the context menu */
            const bool b = mode_->blockSignals(true);
            mode_->setCurrentIndex((int)ppi_->displayMode());
            mode_->blockSignals(b);
        }
        CfarState c;
        {
            Shared &s = shared();
            std::lock_guard<std::mutex> lock(s.mutex);
            c = s.cfar;
        }
        if (!c.known || c.serial == cfarSerial_) return;
        cfarSerial_ = c.serial;
        cfarDefault_dB_ = c.default_dB;
        cfar_->setEnabled(true);
        apply_->setEnabled(true);
        default_->setEnabled(std::isfinite(c.default_dB));
        cfar_->setValue(c.dB);
        cfarInfo_->setText(std::isfinite(c.pfa)
                               ? QStringLiteral("%1: %2 dB (Pfa %3)").arg(detName_).arg(c.dB, 0, 'f', 2).arg(c.pfa, 0, 'e', 1)
                               : QStringLiteral("%1: %2 dB").arg(detName_).arg(c.dB, 0, 'f', 2));
    }

    /* mean time per PRI of every block, from the worker's accumulated times */
    void showTiming()
    {
        ProcessorTiming t;
        uint64_t serial;
        std::string det, trk;
        {
            Shared &s = shared();
            std::lock_guard<std::mutex> lock(s.mutex);
            t = s.timing;
            serial = s.timingSerial;
            det = s.detName;
            trk = s.trackName;
        }
        if (serial == timingSerial_ || t.pris == 0) return;
        timingSerial_ = serial;
        const bool newName = detName_ != QString::fromStdString(det);
        detName_ = QString::fromStdString(det);
        if (newName) cfarSerial_ = 0;                          /* rewrite the CFAR line with the name */
        const double k = 1e3 / (double)t.pris;                 /* ms in total -> us per PRI */
        /* processing thread: the blocks, the console output and the data handed
         * to the PPI; the drawing runs in the GUI thread, at the same time */
        const std::pair<QString, double> rows[] = {
            { QStringLiteral("DBF"), t.dbf_ms }, { QStringLiteral("Matched filter"), t.mf_ms },
            { QStringLiteral("MTI (clutter)"), t.mti_ms }, { QStringLiteral("Doppler FFT+NCI"), t.doppler_ms },
            { detName_, t.cfar_ms }, { QStringLiteral("Plot extractor"), t.extract_ms },
            { QString::fromStdString(trk), t.track_ms }, { QStringLiteral("Console output"), t.console_ms },
            { QStringLiteral("Display (post)"), t.display_ms } };
        QStringList lines;
        double sum = 0.0;
        for (const auto &r : rows) {
            lines << QStringLiteral("%1%2").arg(r.first, -17).arg(r.second * k, 9, 'f', 1);
            sum += r.second * k;
        }
        lines << QStringLiteral("%1%2").arg(QStringLiteral("Total"), -17).arg(sum, 9, 'f', 1);
        lines << QStringLiteral("%1%2").arg(QStringLiteral("PPI drawing (GUI)"), -17).arg(ppi_->guiMs() * k, 9, 'f', 1);
        lines << QStringLiteral("%1 PRIs").arg((qulonglong)t.pris);
        timing_->setText(lines.join(QChar('\n')));
    }

protected:
    void closeEvent(QCloseEvent *event) override
    {
        if (!finished_ && stopFlag_) *stopFlag_ = 1;    /* ask the worker to stop */
        event->accept();
    }

private:
    void requestCfar(double dB)
    {
        if (!std::isfinite(dB)) return;
        {
            Shared &s = shared();
            std::lock_guard<std::mutex> lock(s.mutex);
            s.cfar.requested = true;
            s.cfar.request_dB = dB;
        }
        cfar_->setValue(dB);
        cfarInfo_->setText(QStringLiteral("%1 dB requested ...").arg(dB, 0, 'f', 2));
        ppi_->setFocus();
    }

    PpiWidget             *ppi_;
    volatile sig_atomic_t *stopFlag_;
    QComboBox             *mode_;
    QCheckBox             *tracks_;
    QLabel                *timing_;
    uint64_t               timingSerial_ = 0;
    QString                detName_ = QStringLiteral("CFAR");
    QPushButton           *save_, *apply_, *default_;
    QLabel                *saved_, *cfarInfo_;
    QDoubleSpinBox        *cfar_;
    uint64_t               cfarSerial_ = 0;
    double                 cfarDefault_dB_ = NAN;
    bool                   finished_ = false;
};

/* Without a display (e.g. a Linux server) a snapshot is rendered off screen. */
void use_offscreen_without_display(const char *snapshotPath)
{
#if defined(__linux__)
    if (snapshotPath && qEnvironmentVariableIsEmpty("QT_QPA_PLATFORM") && qEnvironmentVariableIsEmpty("DISPLAY")
        && qEnvironmentVariableIsEmpty("WAYLAND_DISPLAY"))
        qputenv("QT_QPA_PLATFORM", "offscreen");
#else
    (void)snapshotPath;
#endif
}

} /* namespace */

/* ------------------------------------------------------------------------ */
/* Public C interface                                                        */
/* ------------------------------------------------------------------------ */

extern "C" void ppi_post_antenna(double az_deg, double time_s, uint32_t scan)
{
    if (!std::isfinite(az_deg)) return;
    az_deg = wrap360(az_deg);
    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    AntennaState &a = s.antenna;
    if (a.valid && time_s < a.time_s - 1.0) {      /* time went back: a new run */
        a = AntennaState();
        s.newPlots.clear();
        s.clear = true;
    }
    if (a.valid) a.sweep_deg += antenna_step(a, az_deg, scan);
    a.valid = true;
    a.az_deg = az_deg;
    a.time_s = time_s;
    a.scan = scan;
}

extern "C" void ppi_post_plot(const Detection *plot)
{
    if (!plot || !std::isfinite(plot->range_m) || !std::isfinite(plot->az_deg)) return;
    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (s.newPlots.size() < kMaxPending) s.newPlots.push_back(make_mark(*plot, s.antenna, s.nextId++));
}

extern "C" void ppi_post_detections(const Detection *dets, int n, double arcWidth_deg)
{
    if (!dets || n <= 0) return;
    const double arc = std::isfinite(arcWidth_deg) && arcWidth_deg > 0.0 ? std::min(arcWidth_deg, 10.0) : 0.0;
    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    for (int i = 0; i < n && s.newDets.size() < kMaxPendingDets; i++) {
        if (!std::isfinite(dets[i].range_m) || !std::isfinite(dets[i].az_deg)) continue;
        PlotMark m = make_mark(dets[i], s.antenna, s.nextId++);
        m.isDet = true;
        m.arc_deg = arc;
        s.newDets.push_back(m);
    }
}

extern "C" void ppi_post_track(const Track *track)
{
    if (!track) return;
    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.trackUpdates[track->id] = *track;                   /* the GUI needs only the newest state */
}

extern "C" void ppi_post_timing(const ProcessorTiming *timing, const char *detectionName, const char *trackingName)
{
    if (!timing) return;
    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.timing = *timing;
    if (detectionName) s.detName = detectionName;
    if (trackingName) s.trackName = trackingName;
    s.timingSerial++;
}

extern "C" void ppi_post_cfar_threshold(double threshold_dB, double pfa, double default_dB)
{
    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.cfar.known = true;
    s.cfar.dB = threshold_dB;
    s.cfar.pfa = pfa;
    s.cfar.default_dB = default_dB;
    s.cfar.serial++;
}

extern "C" int ppi_take_cfar_request(double *threshold_dB)
{
    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (!s.cfar.requested) return 0;
    s.cfar.requested = false;
    if (threshold_dB) *threshold_dB = s.cfar.request_dB;
    return 1;
}

static std::atomic<int> g_startMode(PPI_SHOW_PLOTS);

extern "C" void ppi_set_display_mode(int mode)
{
    if (mode >= PPI_SHOW_PLOTS && mode <= PPI_SHOW_BOTH) g_startMode.store(mode);
}

extern "C" int ppi_run(int argc, char **argv, double maxRange_m, int (*worker)(void *), void *arg,
                       volatile sig_atomic_t *stopFlag, const char *snapshotPath)
{
    static char defaultName[] = "radar_receiver";
    static char *defaultArgv[] = { defaultName, nullptr };
    int qtArgc = argc;                  /* must outlive the QApplication */
    char **qtArgv = argv;
    if (qtArgc < 1 || !qtArgv || !qtArgv[0]) {
        qtArgc = 1;
        qtArgv = defaultArgv;
    }

    reset_shared();
    use_offscreen_without_display(snapshotPath);
    QApplication app(qtArgc, qtArgv);

    PpiWindow window(maxRange_m, stopFlag);
    PpiWidget &widget = window.ppi();
    window.setWindowTitle(QStringLiteral("Radar PPI"));
    window.resize(1030, 800);
    window.modeBox()->setCurrentIndex(g_startMode.load());   /* DisplayMode has the same order */
    window.show();

    std::atomic<bool> workerDone(false);
    int workerCode = 0;
    std::thread thread([&]() {
        if (worker) workerCode = worker(arg);
        workerDone.store(true);
    });

    bool finishSeen = false;
    QTimer timer;
    QObject::connect(&timer, &QTimer::timeout, [&]() {
        window.refresh();
        if (finishSeen || !workerDone.load()) return;
        finishSeen = true;
        if (snapshotPath) {             /* the snapshot shows all data */
            widget.setPaused(false);
            window.refresh();
        }
        window.setFinished(workerCode);
        if (snapshotPath)
            QCoreApplication::exit(save_snapshot(widget, snapshotPath) ? 0 : 1);
        else if (stopFlag && *stopFlag)
            QCoreApplication::quit();   /* stopped from outside (e.g. Ctrl+C) */
        /* otherwise the window stays open until the user closes it */
    });
    timer.start(kRefreshMs);

    int code = app.exec();
    timer.stop();
    if (!workerDone.load() && stopFlag) *stopFlag = 1;   /* window closed first */
    thread.join();
    if (code == 0) code = workerCode;   /* keep the error code of the worker */
    return code;
}

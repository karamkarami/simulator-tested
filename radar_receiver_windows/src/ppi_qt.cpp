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
 *                          both / raw video), short/long-pulse data, range rings, plot
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
 * plot extractor, CFAR detections before it, both, or the raw video), a button that saves a
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
 * Raw video (display mode Video, ppi_post_video): the DDC output of every PRI
 * averaged over the elements, both receive intervals joined (video.c). The
 * widget keeps the newest data of every 0.25 deg azimuth bin (the mean of the
 * PRIs while the antenna is in the bin) and draws it as an image: colour =
 * power from 3 dB below the noise floor (running 10th percentile of the
 * video) to 50 dB above it.
 * A pixel shows the largest range bin under it. The video stays until the
 * antenna passes again (or Clear).
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
#include <QCursor>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QDateTime>
#include <QHBoxLayout>
#include <QResizeEvent>
#include <QTabWidget>
#include <QImage>
#include <QLabel>
#include <QLinearGradient>
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
constexpr int    kVideoAz     = 1440;     /* video: azimuth bins (0.25 deg)         */
constexpr double kVideoSpan_dB = 50.0;    /* video colours: noise - 3 dB ... + this */
constexpr size_t kMaxPendingVideo = 4000; /* video rows waiting for the GUI thread  */
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

struct VideoRow {                 /* the raw data of one PRI (ppi_post_video, ppi_post_profile) */
    int                kind = 0;      /* 0 video, PROFILE_DBF, PROFILE_MF */
    int                beam = 0;      /* DBF / MF: the beam (0-based)     */
    int                nAvg = 1;      /* complex samples averaged per bin (video: the elements) */
    double             az_deg = 0.0, time_s = 0.0, dR_m = 0.0;
    std::vector<float> v;
};

/* Mean over its 10th percentile [dB] of the mean of n exponential powers (chi-square with
 * 2n degrees of freedom / 2n): x with P(n, x) = 0.1 (regularized lower incomplete gamma,
 * integer n) by bisection, 10 log10(n / x). n = 1: 9.77 dB; 48: 0.86 dB. */
double p10_offset_dB(int n)
{
    if (n < 1) n = 1;
    auto P = [n](double x) {                          /* 1 - e^-x sum_{k<n} x^k / k! */
        double s = 0.0;
        for (int k = 0; k < n; k++) s += std::exp(-x + k * std::log(x) - std::lgamma(k + 1.0));
        return 1.0 - s;
    };
    double lo = 0.0, hi = n;                          /* P(n, n) > 0.1 for every n */
    for (int it = 0; it < 100; it++) {
        const double mid = 0.5 * (lo + hi);
        (P(mid) < 0.1 ? lo : hi) = mid;
    }
    return 10.0 * std::log10(n / (0.5 * (lo + hi)));
}

/* Colours of the video (index 1 ... 255; 0 = no data), as ppi_video_colormap.m. */
QRgb video_colour(int i)
{
    static const double x[7] = { 0, 0.15, 0.35, 0.55, 0.75, 0.9, 1 };
    static const double c[7][3] = { { 0, 0, 0 }, { 0, 0, 0.55 }, { 0, 0.6, 1 }, { 0, 1, 0.3 },
                                     { 1, 1, 0 }, { 1, 0.25, 0 }, { 1, 1, 1 } };
    const double t = (i - 1) / 254.0;
    int k = 0;
    while (k < 5 && t > x[k + 1]) k++;
    const double f = (t - x[k]) / (x[k + 1] - x[k]);
    auto ch = [&](int j) { return (int)std::lround(255.0 * (c[k][j] + f * (c[k + 1][j] - c[k][j]))); };
    return qRgb(ch(0), ch(1), ch(2));
}

std::atomic<int> g_videoWanted(0);   /* the PPI shows the video: the worker computes it */
std::atomic<int> g_profileKind(0);   /* PROFILE_DBF / PROFILE_MF while the PPI shows that output */
std::atomic<int> g_beam(0);          /* the beam of the DBF / MF output and the range-Doppler map */
std::atomic<int> g_rdWanted(0);      /* the range-Doppler tab is shown */
std::vector<double> &beam_elevations()   /* the beams of the beam selector (ppi_set_beams) */
{
    static std::vector<double> el(1, NAN);
    return el;
}

/* Range-Doppler map of one beam (ppi_post_rdmap): the copy for the GUI thread. */
struct RdFrame {
    bool                valid = false;
    int                 beam = 0, nBins = 0, nfft = 0, notchBins = 0;
    double              dR_m = 0.0, noise = 0.0, time_s = 0.0, az_deg = 0.0;
    std::vector<float>  power;        /* [nfft][nBins] */
    std::vector<double> vel;          /* [nfft]        */
    std::vector<Detection> dets;
};

struct Shared {
    std::mutex            mutex;
    AntennaState          antenna;
    std::vector<PlotMark> newPlots;   /* plots not yet taken by the GUI thread */
    std::vector<PlotMark> newDets;    /* CFAR detections not yet taken         */
    std::vector<VideoRow> newVideo;   /* video rows not yet taken              */
    RdFrame               rd;         /* newest range-Doppler map              */
    uint64_t              rdSerial = 0;
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
    s.newVideo.clear();
    s.rd = RdFrame();
    s.rdSerial = 0;
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

/* A row of the timing panel: the label shortened to 15 characters ("1-D" ->
 * "1D", then "range" -> "rng", "Doppler" -> "Dop", then cut), the value
 * right-aligned in 9 (24 characters in all). */
QString timing_label(QString name)
{
    name.replace(QStringLiteral("-D"), QStringLiteral("D"));
    if (name.length() > 15)
        name.replace(QStringLiteral(" range"), QStringLiteral(" rng")).replace(QStringLiteral(" Doppler"), QStringLiteral(" Dop"));
    return name.length() > 15 ? name.left(14) + QChar('.') : name;
}

QString timing_row(const QString &name, double us)
{
    return QStringLiteral("%1%2").arg(timing_label(name), -15).arg(us, 9, 'f', 1);
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
/* after / before the plot extractor, both, raw video, DBF output and pulse compression
 * output of the selected beam (the order of the panel's box and of PPI_SHOW_...) */
enum class DisplayMode { Plots, Detections, Both, Video, DbfOut, MfOut };

/* The kind of raw data a display mode shows: 0 video, PROFILE_DBF, PROFILE_MF; -1 none. */
int raw_kind(DisplayMode m)
{
    return m == DisplayMode::Video ? 0 : m == DisplayMode::DbfOut ? PROFILE_DBF : m == DisplayMode::MfOut ? PROFILE_MF : -1;
}

QString beam_text(int beam)
{
    const std::vector<double> &el = beam_elevations();
    const double e = beam >= 0 && beam < (int)el.size() ? el[beam] : NAN;
    return std::isfinite(e) ? QStringLiteral("beam %1 (el %2)").arg(beam + 1).arg(deg_text(e, 1))
                            : QStringLiteral("beam %1").arg(beam + 1);
}

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
    void   setDisplayMode(DisplayMode m)
    {
        if (raw_kind(m) != raw_kind(mode_)) clearVideo();  /* other raw data */
        mode_ = m;
        g_videoWanted.store(m == DisplayMode::Video);
        g_profileKind.store(raw_kind(m) > 0 ? raw_kind(m) : PROFILE_NONE);
        videoImageDirty_ = true;
        update();
    }
    /* beam (0-based) of the DBF / MF output and of the range-Doppler map */
    void   setBeam(int beam)
    {
        if (beam == g_beam.load()) return;
        g_beam.store(beam);
        if (raw_kind(mode_) > 0) clearVideo();
        update();
    }
    int    beam() const { return g_beam.load(); }
    /* video: azimuth bins with data, noise floor [dB], colour index of a ground point (0 = none) */
    int    videoBinsFilled() const
    {
        int n = 0;
        for (int a = 0; a < kVideoAz && videoBins_ > 0; a++)
            n += *std::max_element(videoIdx_.begin() + (size_t)a * videoBins_,
                                   videoIdx_.begin() + (size_t)(a + 1) * videoBins_) > 0;
        return n;
    }
    double videoNoise_dB() const { return videoNoise_dB_; }
    int    videoIndexAt(double east_m, double north_m) const
    {
        if (videoBins_ < 1) return 0;
        const int b = (int)std::floor(std::hypot(east_m, north_m) / videoDR_);
        const int a = std::min(kVideoAz - 1, (int)std::floor(wrap360(std::atan2(east_m, north_m) / kDeg) / (360.0 / kVideoAz)));
        return b >= 0 && b < videoBins_ ? videoIdx_[(size_t)a * videoBins_ + b] : 0;
    }
    DisplayMode displayMode() const { return mode_; }
    QString cursorText() const { return cursorText_; }
    void   setShowShort(bool on) { showShort_ = on; update(); }
    void   setShowLong(bool on) { showLong_ = on; update(); }
    void   setMinSnr(double dB) { minSnr_dB_ = dB; update(); }
    void   setHistory(History h) { history_ = h; refresh(); }
    void   setPaused(bool on) { paused_ = on; update(); }
    bool   paused() const { return paused_; }
    void   clearPlots() { plots_.clear(); dets_.clear(); tracks_.clear(); pinned_.clear(); clearVideo(); update(); }
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
    void   drawVideo(QPainter &p);
    void   addVideoRow(const VideoRow &row);
    void   quantiseVideoRow(int az);
    void   clearVideo();
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
    /* raw video */
    std::vector<VideoRow> incomingVideo_;
    int                   videoBins_ = 0;
    double                videoDR_ = 0.0;
    std::vector<float>    videoMap_;      /* mean power [kVideoAz][videoBins_] (0 = no data) */
    std::vector<uint8_t>  videoIdx_;      /* its colour index (0 = no data)            */
    std::vector<float>    videoSum_;      /* PRIs of the current azimuth bin            */
    std::vector<float>    videoTmp_;
    int                   videoN_ = 0, videoAzBin_ = -1;
    double                videoNoise_dB_ = NAN, videoLo_dB_ = NAN;   /* colour scale from videoLo_dB_ */
    QImage                videoImage_;
    bool                  videoImageDirty_ = true;
    std::vector<int32_t>  lutCell_;       /* per pixel: az * bins + first range bin, -1 outside */
    int                   lutK_ = 1;      /* range bins per pixel                       */
    QSize                 lutSize_;
    double                lutZoom_ = 0.0, lutE_ = 0.0, lutN_ = 0.0;
    int                   lutBins_ = 0;
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
        incomingVideo_.swap(s.newVideo);
        s.newVideo.clear();
        antenna_ = s.antenna;
        clear = s.clear;
        s.clear = false;
    }
    if (clear) { plots_.clear(); dets_.clear(); tracks_.clear(); pinned_.clear(); clearVideo(); }
    for (const VideoRow &row : incomingVideo_)           /* only the raw data of this mode (and beam) */
        if (row.kind == raw_kind(mode_) && (row.kind == 0 || row.beam == g_beam.load())) addVideoRow(row);
    incomingVideo_.clear();
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
    if (raw_kind(mode_) >= 0) return false;            /* raw data instead of plots and detections */
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
                                       : mode_ == DisplayMode::Detections ? DisplayMode::Both
                                       : mode_ == DisplayMode::Both ? DisplayMode::Video
                                       : mode_ == DisplayMode::Video ? DisplayMode::DbfOut
                                       : mode_ == DisplayMode::DbfOut ? DisplayMode::MfOut : DisplayMode::Plots);
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
        { "Plots and detections", DisplayMode::Both },
        { "Video (DDC output, elements averaged)", DisplayMode::Video },
        { "DBF out (selected beam)", DisplayMode::DbfOut },
        { "MF out (pulse compression, selected beam)", DisplayMode::MfOut } };
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
    rows << QStringLiteral("<b>Track %1 (%2, %3-D %4)</b>").arg(t.id).arg(t.confirmed ? "confirmed" : "tentative").arg(t.dim)
                .arg(t.method == TRACK_METHOD_PMBM ? "PMBM" : "EKF")
         << QStringLiteral("Range: %1 km").arg(R / 1e3, 0, 'f', 3)
         << QStringLiteral("Azimuth: %1").arg(deg_text(wrap360(std::atan2(t.pos_m[0], t.pos_m[1]) / kDeg), 2));
    if (t.dim == 3)
        rows << QStringLiteral("Height: %1 m").arg(t.alt_m, 0, 'f', 0)
             << QStringLiteral("Climb: %1 m/s").arg(t.vel_mps[2], 0, 'f', 1);
    rows << QStringLiteral("Speed: %1 m/s").arg(std::hypot(t.vel_mps[0], t.vel_mps[1]), 0, 'f', 1)
         << QStringLiteral("Heading: %1").arg(deg_text(wrap360(std::atan2(t.vel_mps[0], t.vel_mps[1]) / kDeg), 1));
    const double acc = std::hypot(t.acc_mps2[0], t.acc_mps2[1]);
    if (acc > 0.0) rows << QStringLiteral("Acceleration: %1 m/s\u00b2").arg(acc, 0, 'f', 2);   /* CA model */
    if (t.method == TRACK_METHOD_PMBM) {                     /* existence, IMM model probabilities */
        static const char *names[3] = { "CV", "CA", "CT" };
        QStringList modes;
        for (int j = 0; j < t.nModes; j++)
            modes << QStringLiteral("%1 %2").arg(names[t.modes[j]]).arg(t.modeProb[j], 0, 'f', 2);
        rows << QStringLiteral("Existence: %1").arg(t.existence, 0, 'f', 3)
             << QStringLiteral("Models: %1").arg(modes.join(QStringLiteral(" ")));
    }
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

/* ---- raw video ---------------------------------------------------------- */

void PpiWidget::clearVideo()
{
    std::fill(videoMap_.begin(), videoMap_.end(), 0.0f);
    std::fill(videoIdx_.begin(), videoIdx_.end(), (uint8_t)0);
    videoN_ = 0;
    videoAzBin_ = -1;
    videoImageDirty_ = true;
}

/* Colour indices of azimuth bin az from its mean power (scale from videoLo_dB_). */
void PpiWidget::quantiseVideoRow(int az)
{
    const float *m = &videoMap_[(size_t)az * videoBins_];
    uint8_t *q = &videoIdx_[(size_t)az * videoBins_];
    const double k = 254.0 / (kVideoSpan_dB + 3.0);
    for (int b = 0; b < videoBins_; b++) {
        if (!(m[b] > 0.0f)) { q[b] = 0; continue; }
        const double i = 1.0 + std::round((10.0 * std::log10((double)m[b]) - videoLo_dB_) * k);
        q[b] = (uint8_t)std::max(1.0, std::min(255.0, i));
    }
}

/* The video of one PRI into the bin of its azimuth (mean of the PRIs while the
 * antenna is in the bin; bins skipped by a jump below 10 deg get the same data). */
void PpiWidget::addVideoRow(const VideoRow &row)
{
    const int nB = (int)row.v.size();
    if (nB < 1 || !(row.dR_m > 0.0)) return;
    if (nB != videoBins_ || row.dR_m != videoDR_) {       /* first row, or a new layout */
        videoBins_ = nB;
        videoDR_ = row.dR_m;
        videoMap_.assign((size_t)kVideoAz * nB, 0.0f);
        videoIdx_.assign((size_t)kVideoAz * nB, 0);
        videoSum_.assign(nB, 0.0f);
        videoN_ = 0;
        videoAzBin_ = -1;
        videoNoise_dB_ = videoLo_dB_ = NAN;
        lutBins_ = 0;                                     /* rebuild the pixel lookup */
    }
    /* noise floor: running 10th percentile of the video in dB (the raw video is
     * mostly clutter: the long pulse smears near clutter over its length) */
    videoTmp_.clear();
    for (float x : row.v) if (x > 0.0f) videoTmp_.push_back(x);
    if (!videoTmp_.empty()) {
        auto mid = videoTmp_.begin() + (videoTmp_.size() - 1) / 10;
        std::nth_element(videoTmp_.begin(), mid, videoTmp_.end());
        /* the 10th percentile of a bin is below the noise power: DBF / MF out one complex
         * sample per bin (exponential power: 9.77 dB), the video the mean over the elements
         * (one element: also 9.77 dB, 48 elements: 0.86 dB) */
        const double m = 10.0 * std::log10((double)*mid) + p10_offset_dB(row.kind == 0 ? row.nAvg : 1);
        videoNoise_dB_ = std::isfinite(videoNoise_dB_) ? 0.98 * videoNoise_dB_ + 0.02 * m : m;
    }
    bool requantise = false;
    if (std::isfinite(videoNoise_dB_) && !(std::fabs(videoNoise_dB_ - 3.0 - videoLo_dB_) < 0.5)) {
        videoLo_dB_ = videoNoise_dB_ - 3.0;               /* the scale moved: all bins again */
        requantise = true;
    }
    const int b = std::min(kVideoAz - 1, (int)std::floor(wrap360(row.az_deg) / (360.0 / kVideoAz)));
    if (b == videoAzBin_) {
        for (int i = 0; i < nB; i++) videoSum_[i] += row.v[i];
        videoN_++;
    } else {
        if (videoAzBin_ >= 0) {
            const int gap = (b - videoAzBin_ + kVideoAz) % kVideoAz - 1;
            for (int g = 1; g <= gap && gap * (360.0 / kVideoAz) < 10.0; g++) {
                const int a = (videoAzBin_ + g) % kVideoAz;
                std::copy(row.v.begin(), row.v.end(), videoMap_.begin() + (size_t)a * nB);
                if (!requantise) quantiseVideoRow(a);
            }
        }
        videoAzBin_ = b;
        std::copy(row.v.begin(), row.v.end(), videoSum_.begin());
        videoN_ = 1;
    }
    float *m = &videoMap_[(size_t)b * nB];
    for (int i = 0; i < nB; i++) m[i] = videoSum_[i] / (float)videoN_;
    if (requantise)
        for (int a = 0; a < kVideoAz; a++) quantiseVideoRow(a);
    else
        quantiseVideoRow(b);
    videoImageDirty_ = true;
}

/* The video as an image of the widget: per pixel the largest colour index of
 * the range bins under it (a lookup table per pixel, rebuilt when the view
 * changes; the image is rebuilt when new video arrived). */
void PpiWidget::drawVideo(QPainter &p)
{
    if (raw_kind(mode_) < 0 || videoBins_ < 1) return;
    const QSize sz = size();
    if (sz.isEmpty()) return;
    if (lutSize_ != sz || lutZoom_ != zoom_ || lutE_ != viewE_ || lutN_ != viewN_ || lutBins_ != videoBins_) {
        lutSize_ = sz;  lutZoom_ = zoom_;  lutE_ = viewE_;  lutN_ = viewN_;  lutBins_ = videoBins_;
        const double pixel_m = 1.0 / scale();
        lutK_ = std::max(1, (int)std::ceil(pixel_m / videoDR_));
        lutCell_.assign((size_t)sz.width() * sz.height(), -1);
        const double rMax = std::min(maxRange_m_, videoBins_ * videoDR_);
        for (int y = (int)viewRect_.top(); y < sz.height(); y++) {
            for (int x = 0; x < sz.width(); x++) {
                const QPointF g = toGround(QPointF(x + 0.5, y + 0.5));
                const double r = std::hypot(g.x(), g.y());
                if (r >= rMax) continue;
                const int b0 = std::max(0, std::min(videoBins_ - lutK_, (int)std::floor(r / videoDR_ - 0.5 * (lutK_ - 1))));
                const int a = std::min(kVideoAz - 1, (int)std::floor(wrap360(std::atan2(g.x(), g.y()) / kDeg) / (360.0 / kVideoAz)));
                lutCell_[(size_t)y * sz.width() + x] = a * videoBins_ + std::max(0, b0);
            }
        }
        videoImageDirty_ = true;
    }
    if (videoImageDirty_ || videoImage_.size() != sz) {
        static QRgb colours[256];
        static bool init = false;
        if (!init) {
            colours[0] = qRgba(0, 0, 0, 0);
            for (int i = 1; i < 256; i++) colours[i] = video_colour(i);
            init = true;
        }
        videoImage_ = QImage(sz, QImage::Format_ARGB32_Premultiplied);
        const int K = std::min(lutK_, videoBins_);
        const uint8_t *q = videoIdx_.data();
        for (int y = 0; y < sz.height(); y++) {
            QRgb *line = reinterpret_cast<QRgb *>(videoImage_.scanLine(y));
            const int32_t *cell = &lutCell_[(size_t)y * sz.width()];
            for (int x = 0; x < sz.width(); x++) {
                uint8_t v = 0;
                if (cell[x] >= 0)
                    for (int k = 0; k < K; k++) v = std::max(v, q[cell[x] + k]);
                line[x] = colours[v];
            }
        }
        videoImageDirty_ = false;
    }
    p.drawImage(QPointF(0.0, 0.0), videoImage_);
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
        if (raw_kind(mode_) >= 0) {
            const QString what = mode_ == DisplayMode::Video ? QStringLiteral("Video")
                               : QStringLiteral("%1 %2").arg(mode_ == DisplayMode::DbfOut ? "DBF out" : "MF out").arg(beam_text(g_beam.load()));
            line1 += std::isfinite(videoNoise_dB_) ? QStringLiteral("    %1 (noise %2 dB)").arg(what).arg(videoNoise_dB_, 0, 'f', 1)
                                                   : QStringLiteral("    %1 (no data yet)").arg(what);
        } else if (mode_ != DisplayMode::Detections)
            line1 += QStringLiteral("    Plots %1 / %2").arg((qulonglong)visibleCount()).arg((qulonglong)plots_.size());
        if (mode_ != DisplayMode::Plots && raw_kind(mode_) < 0)
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
    enum Kind { PlotS, PlotL, DetS, DetL, VideoK, TrackK, Antenna, Help };
    std::vector<std::pair<Kind, QString>> rows;
    if (raw_kind(mode_) >= 0) {
        rows.emplace_back(Help, mode_ == DisplayMode::Video ? QStringLiteral("Video [dB above noise], elements averaged")
                                : mode_ == DisplayMode::DbfOut ? QStringLiteral("DBF out [dB above noise], %1").arg(beam_text(g_beam.load()))
                                : QStringLiteral("MF out [dB above noise], %1").arg(beam_text(g_beam.load())));
        rows.emplace_back(VideoK, QString());
        rows.emplace_back(Help, QString());              /* (room for the scale labels) */
    } else if (mode_ != DisplayMode::Detections) {
        rows.emplace_back(PlotS, QStringLiteral("Plot, short pulse"));
        rows.emplace_back(PlotL, QStringLiteral("Plot, long pulse"));
    }
    if (mode_ != DisplayMode::Plots && raw_kind(mode_) < 0) {
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
    if (raw_kind(mode_) >= 0) textWidth = std::max(textWidth, 230);
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
        } else if (kind == VideoK) {                       /* colour scale of the video */
            const double w = 26.0 + textWidth - 48.0, x0 = box.left() + pad + 8.0;
            QLinearGradient g(x0, 0.0, x0 + w, 0.0);
            for (int k = 1; k <= 255; k += 6) g.setColorAt((k - 1) / 254.0, QColor(video_colour(k)));
            g.setColorAt(1.0, QColor(video_colour(255)));
            p.setPen(Qt::NoPen);
            p.setBrush(g);
            p.drawRect(QRectF(x0, y - 7.0, w, 7.0));
            p.setPen(QColor(Qt::white));
            for (int dB = 0; dB <= (int)kVideoSpan_dB; dB += 10) {
                const double x = x0 + w * (dB + 3.0) / (kVideoSpan_dB + 3.0);
                const QString t = QString::number(dB);
                p.drawText(QPointF(x - fm.horizontalAdvance(t) / 2.0, y + fm.ascent() + 1.0), t);
            }
            continue;
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
    drawVideo(p);
    drawAntenna(p);
    drawDetections(p);
    drawPlots(p);
    drawTracks(p);
    p.restore();
    drawPins(p);
    drawStatus(p);
    drawLegend(p);
}


/* ------------------------------------------------------------------------ */
/* Range-Doppler map of the selected beam (second tab)                       */
/* ------------------------------------------------------------------------ */

/* "Nice" axis step (1, 2, 5 x 10^k) for about n ticks over a span. */
double nice_step(double span, int n)
{
    const double raw = span / std::max(1, n), p = std::pow(10.0, std::floor(std::log10(raw)));
    const double f = raw / p;
    return (f < 1.5 ? 1.0 : f < 3.5 ? 2.0 : f < 7.5 ? 5.0 : 10.0) * p;
}

/* What the plot below the map shows: the Doppler spectrum of one range bin or
 * the range profile of one Doppler bin (the panel's box; --ppi-rd-plot). */
enum class RdLower { Spectrum, RangeProfile };

/* The power map of the beam (range x radial velocity, colour = dB above the
 * noise floor of the map, from -3 dB to the strongest cell), the CFAR
 * detections of the beam (white circles), the MTI notch (dashed), and below
 * it a cut through the map: the Doppler spectrum of one range bin or the range
 * profile of one Doppler bin, through the cell under the cursor, else the
 * selected cell (a click), else the strongest detection (or strongest cell).
 * Both plots zoom (wheel at the cursor; Shift: only the horizontal axis,
 * Ctrl: only the vertical one; keys + -), pan (drag; arrow keys), reset
 * (double click; key R), show a data tip of the cell, detection or point under
 * the cursor and pin it with a click (click it again, or key C: remove). The
 * view and the data tips stay while new maps arrive. */
class RdWidget : public QWidget {
public:
    RdWidget()
    {
        setMouseTracking(true);
        setMinimumSize(300, 250);
        setAttribute(Qt::WA_OpaquePaintEvent);
        setFocusPolicy(Qt::ClickFocus);
    }
    void setFrame(const RdFrame &f)
    {
        const bool sameLayout = frame_.valid && f.valid && f.nBins == frame_.nBins && f.nfft == frame_.nfft &&
                                f.dR_m == frame_.dR_m && f.beam == frame_.beam;
        frame_ = f;
        image_ = QImage();
        if (!sameLayout) { mapAuto_ = true; lowAuto_ = true; pins_.clear(); selB_ = selD_ = -1; }
        if (!frame_.valid) { first_ = 0; last_ = -1; update(); return; }
        /* range bins with data, dB above noise, the scale */
        first_ = frame_.nBins; last_ = -1;
        maxDb_ = -1e9;
        const double noise = frame_.noise > 0 ? frame_.noise : 1.0;
        db_.assign(frame_.power.size(), NAN);
        for (int d = 0; d < frame_.nfft; d++)
            for (int b = 0; b < frame_.nBins; b++) {
                const float v = frame_.power[(size_t)d * frame_.nBins + b];
                if (!(v > 0)) continue;
                const double x = 10.0 * std::log10(v / noise);
                db_[(size_t)d * frame_.nBins + b] = (float)x;
                first_ = std::min(first_, b);
                last_ = std::max(last_, b);
                maxDb_ = std::max(maxDb_, x);
            }
        hiDb_ = std::max(20.0, std::ceil(maxDb_ / 5.0) * 5.0);
        /* the cut without a cursor or a selection: strongest detection, else strongest cell */
        autoB_ = autoD_ = -1;
        double best = -1e9;
        for (const Detection &dt : frame_.dets) {
            const int b = binOf(dt.range_m);
            if (b >= 0 && b < frame_.nBins && dt.dIdx >= 0 && dt.dIdx < frame_.nfft && dt.power > best) {
                best = dt.power; autoB_ = b; autoD_ = dt.dIdx;
            }
        }
        if (autoB_ < 0) {
            for (size_t i = 0; i < db_.size(); i++)
                if (db_[i] == db_[i] && db_[i] > best) {
                    best = db_[i]; autoB_ = (int)(i % frame_.nBins); autoD_ = (int)(i / frame_.nBins);
                }
        }
        if (hasFrame()) fitViews();
        update();
    }
    const RdFrame &frame() const { return frame_; }
    bool   hasFrame() const { return frame_.valid && last_ >= first_; }
    double maxDb() const { return maxDb_; }
    /* the cut of the lower plot: range bin of the spectrum, Doppler bin of the range profile */
    int    spectrumBin() const { return curB_ >= 0 ? curB_ : selB_ >= 0 ? selB_ : autoB_; }
    int    profileBin() const { return curD_ >= 0 ? curD_ : selD_ >= 0 ? selD_ : autoD_; }
    RdLower lower() const { return lower_; }
    void   setLower(RdLower m) { if (m != lower_) { lower_ = m; lowAuto_ = true; fitViews(); update(); } }
    void   setPaused(bool on) { if (on != paused_) { paused_ = on; update(); } }
    QRectF mapRect() const { return mapRect_; }
    QRectF lowerRect() const { return lowRect_; }
    QString readout() const { return readout_; }
    QString tipText() const { return tip_; }            /* the data tip under the cursor (tests) */
    int    pinCount() const { return (int)pins_.size(); }
    /* the visible part of the map [range bin, Doppler bin coordinates] and of the lower plot */
    QRectF mapView() const { return QRectF(QPointF(vx0_, vy0_), QPointF(vx1_, vy1_)); }
    QRectF lowerView() const { return QRectF(QPointF(lx0_, ly0_), QPointF(lx1_, ly1_)); }
    /* screen point of a range [m] and Doppler bin (the centre of its row) */
    QPointF screenOf(double range_m, int dBin) const
    {
        return mapScreen(range_m / frame_.dR_m, dBin + 0.5);
    }
    /* screen point of the lower plot at a range bin (range profile) or Doppler bin (spectrum) and a level [dB] */
    QPointF lowerScreenOf(int bin, double db) const { return lowScreen(lowerX(bin), db); }
    void hoverAt(QPointF pos) { layoutRects(); updateCursor(pos); }
    void resetViews() { mapAuto_ = lowAuto_ = true; fitViews(); update(); }

protected:
    void mouseMoveEvent(QMouseEvent *e) override;
    void mousePressEvent(QMouseEvent *e) override;
    void mouseReleaseEvent(QMouseEvent *e) override;
    void mouseDoubleClickEvent(QMouseEvent *e) override;
    void wheelEvent(QWheelEvent *e) override;
    void keyPressEvent(QKeyEvent *e) override;
    void contextMenuEvent(QContextMenuEvent *e) override;
    void leaveEvent(QEvent *) override
    {
        curB_ = curD_ = -1; readout_.clear(); tip_.clear(); tipKey_ = -2; QToolTip::hideText(); update();
    }
    void resizeEvent(QResizeEvent *) override { image_ = QImage(); update(); }
    void paintEvent(QPaintEvent *) override;

private:
    struct Pin { int b, d; };                          /* a pinned cell (range bin, Doppler bin) */
    enum Area { None, Map, Lower };

    int binOf(double range_m) const { return (int)std::floor(range_m / frame_.dR_m); }
    double dv() const { return frame_.nfft > 1 ? frame_.vel[1] - frame_.vel[0] : 1.0; }
    /* velocity of a Doppler row coordinate (row d spans [d, d+1)) and back */
    double velAt(double y) const { return frame_.vel[0] + (y - 0.5) * dv(); }
    double rowOf(double v) const { return (v - frame_.vel[0]) / dv() + 0.5; }
    double cellDb(int b, int d) const
    {
        return b >= 0 && b < frame_.nBins && d >= 0 && d < frame_.nfft ? db_[(size_t)d * frame_.nBins + b] : NAN;
    }
    /* lower plot: x of a bin of the cut (velocity [m/s] of a Doppler bin, or range [km] of a range bin) */
    double lowerX(int bin) const
    {
        return lower_ == RdLower::Spectrum ? frame_.vel[std::max(0, std::min(frame_.nfft - 1, bin))]
                                           : (bin + 0.5) * frame_.dR_m / 1e3;
    }
    void lowerFull(double &x0, double &x1, double &y0, double &y1) const
    {
        if (lower_ == RdLower::Spectrum) {
            const double a = frame_.vel[0], b = frame_.vel[frame_.nfft - 1], h = 0.5 * std::fabs(dv());
            x0 = std::min(a, b) - h;  x1 = std::max(a, b) + h;
        } else {
            x0 = first_ * frame_.dR_m / 1e3;  x1 = (last_ + 1) * frame_.dR_m / 1e3;
        }
        y0 = -10.0;  y1 = std::max(10.0, std::ceil((maxDb_ + 3.0) / 5.0) * 5.0);   /* room above the peak */
    }
    void fitViews()
    {
        if (!hasFrame()) return;
        const double fx0 = first_, fx1 = last_ + 1.0, fy0 = 0.0, fy1 = frame_.nfft;
        if (mapAuto_) { vx0_ = fx0; vx1_ = fx1; vy0_ = fy0; vy1_ = fy1; }
        else clampMap();
        double a, b, c, d;
        lowerFull(a, b, c, d);
        if (lowAuto_) { lx0_ = a; lx1_ = b; ly0_ = c; ly1_ = d; }
        else clampLower();
    }
    static void clampSpan(double &lo, double &hi, double fullLo, double fullHi, double minSpan)
    {
        double span = std::max(minSpan, std::min(hi - lo, fullHi - fullLo));
        const double mid = 0.5 * (lo + hi);
        lo = mid - 0.5 * span;  hi = mid + 0.5 * span;
        if (lo < fullLo) { hi += fullLo - lo; lo = fullLo; }
        if (hi > fullHi) { lo -= hi - fullHi; hi = fullHi; }
    }
    void clampMap()
    {
        clampSpan(vx0_, vx1_, first_, last_ + 1.0, std::min(4.0, last_ + 1.0 - first_));
        clampSpan(vy0_, vy1_, 0.0, frame_.nfft, std::min(2.0, (double)frame_.nfft));
    }
    void clampLower()
    {
        double a, b, c, d;
        lowerFull(a, b, c, d);
        const double minX = lower_ == RdLower::Spectrum ? 2.0 * std::fabs(dv()) : 4.0 * frame_.dR_m / 1e3;
        clampSpan(lx0_, lx1_, a, b, std::min(minX, b - a));
        clampSpan(ly0_, ly1_, c - 30.0, d + 30.0, 3.0);
    }
    void layoutRects()
    {
        const double w = width(), h = height(), left = 64.0, right = 78.0, top = 50.0;
        const double gap = 58.0, bottom = 40.0, lowH = std::max(80.0, 0.28 * h);
        mapRect_ = QRectF(left, top, std::max(50.0, w - left - right), std::max(50.0, h - top - gap - lowH - bottom));
        lowRect_ = QRectF(left, mapRect_.bottom() + gap, mapRect_.width(), lowH);
    }
    QPointF mapScreen(double x, double y) const
    {
        return QPointF(mapRect_.left() + (x - vx0_) / (vx1_ - vx0_) * mapRect_.width(),
                       mapRect_.top() + (y - vy0_) / (vy1_ - vy0_) * mapRect_.height());
    }
    QPointF mapData(QPointF s) const
    {
        return QPointF(vx0_ + (s.x() - mapRect_.left()) / mapRect_.width() * (vx1_ - vx0_),
                       vy0_ + (s.y() - mapRect_.top()) / mapRect_.height() * (vy1_ - vy0_));
    }
    QPointF lowScreen(double x, double y) const
    {
        return QPointF(lowRect_.left() + (x - lx0_) / (lx1_ - lx0_) * lowRect_.width(),
                       lowRect_.bottom() - (y - ly0_) / (ly1_ - ly0_) * lowRect_.height());
    }
    QPointF lowData(QPointF s) const
    {
        return QPointF(lx0_ + (s.x() - lowRect_.left()) / lowRect_.width() * (lx1_ - lx0_),
                       ly0_ + (lowRect_.bottom() - s.y()) / lowRect_.height() * (ly1_ - ly0_));
    }
    Area areaAt(QPointF pos) const
    {
        if (!hasFrame()) return None;
        return mapRect_.contains(pos) ? Map : lowRect_.contains(pos) ? Lower : None;
    }
    /* the points of the lower plot: (x, dB, bin) of the cut, sorted by x */
    struct LowPoint { double x, db; int bin; };
    std::vector<LowPoint> lowerPoints() const
    {
        std::vector<LowPoint> pts;
        if (lower_ == RdLower::Spectrum) {
            const int sb = spectrumBin();
            if (sb < first_ || sb > last_) return pts;
            for (int d = 0; d < frame_.nfft; d++) {
                const double v = cellDb(sb, d);
                if (v == v) pts.push_back({ lowerX(d), v, d });
            }
        } else {
            const int pd = profileBin();
            if (pd < 0 || pd >= frame_.nfft) return pts;
            for (int b = first_; b <= last_; b++) {
                const double v = cellDb(b, pd);
                if (v == v) pts.push_back({ lowerX(b), v, b });
            }
        }
        std::sort(pts.begin(), pts.end(), [](const LowPoint &a, const LowPoint &b) { return a.x < b.x; });
        return pts;
    }
    /* the cell of a bin of the lower plot */
    Pin lowerCell(int bin) const
    {
        return lower_ == RdLower::Spectrum ? Pin{ spectrumBin(), bin } : Pin{ bin, profileBin() };
    }
    const Detection *detectionAt(QPointF pos) const
    {
        const Detection *best = nullptr;
        double bestD = 8.0;
        for (const Detection &dt : frame_.dets) {
            const QPointF d = screenOf(dt.range_m, dt.dIdx) - pos;
            const double dist = std::hypot(d.x(), d.y());
            if (dist <= bestD) { bestD = dist; best = &dt; }
        }
        return best;
    }
    const Detection *detectionOf(int b, int d) const
    {
        for (const Detection &dt : frame_.dets)
            if (binOf(dt.range_m) == b && dt.dIdx == d) return &dt;
        return nullptr;
    }
    /* data tip of a cell (with its detection, if any) */
    QString cellText(int b, int d, bool html) const
    {
        const double v = cellDb(b, d);
        QStringList rows;
        const Detection *dt = detectionOf(b, d);
        rows << (dt ? QStringLiteral("CFAR detection, %1 pulse").arg(dt->pulse == PULSE_LONG ? "long" : "short")
                    : QStringLiteral("Range-Doppler cell"))
             << QStringLiteral("Range: %1 km").arg((b + 0.5) * frame_.dR_m / 1e3, 0, 'f', 3)
             << QStringLiteral("Radial velocity: %1 m/s").arg(frame_.vel[d], 0, 'f', 2)
             << (v == v ? QStringLiteral("Power: %1 dB above noise").arg(v, 0, 'f', 1) : QStringLiteral("Power: no cell"));
        if (dt) rows << QStringLiteral("SNR: %1 dB").arg(dt->snr_dB, 0, 'f', 1);
        rows << QStringLiteral("Range bin %1, Doppler bin %2").arg(b).arg(d);
        if (html) rows[0] = QStringLiteral("<b>%1</b>").arg(rows[0]);
        return rows.join(html ? QStringLiteral("<br>") : QStringLiteral("\n"));
    }
    void updateCursor(QPointF pos);
    void zoomAt(Area a, QPointF pos, double fx, double fy);
    void panBy(Area a, QPointF pixels);
    void togglePin(QPointF pos);
    void drawPins(QPainter &p) const;
    void drawLower(QPainter &p, const QFontMetrics &fm);
    static QRgb colour_db(double x, double lo, double hi)
    {
        const int i = 1 + (int)std::lround(std::max(0.0, std::min(1.0, (x - lo) / (hi - lo))) * 254.0);
        return video_colour(i);
    }

    RdFrame            frame_;
    std::vector<float> db_;
    int                first_ = 0, last_ = -1, autoB_ = -1, autoD_ = -1;
    int                curB_ = -1, curD_ = -1;       /* the cell under the cursor (map) */
    int                selB_ = -1, selD_ = -1;       /* the selected cell (click on the map) */
    double             maxDb_ = 0.0, hiDb_ = 20.0;
    double             vx0_ = 0, vx1_ = 1, vy0_ = 0, vy1_ = 1;     /* map view: range bins, Doppler rows */
    double             lx0_ = 0, lx1_ = 1, ly0_ = 0, ly1_ = 1;     /* lower plot view: x, dB */
    bool               mapAuto_ = true, lowAuto_ = true, paused_ = false;
    RdLower            lower_ = RdLower::Spectrum;
    std::vector<Pin>   pins_;
    QImage             image_;
    QRectF             mapRect_, lowRect_;
    QString            readout_, tip_;
    long long          tipKey_ = -2;                 /* what the shown data tip belongs to */
    Area               dragArea_ = None;
    QPointF            pressPos_, lastPos_;
    bool               dragging_ = false;
};

void RdWidget::updateCursor(QPointF pos)
{
    curB_ = curD_ = -1;
    readout_.clear();
    QString tip;
    long long key = -2;
    const Area a = areaAt(pos);
    if (a == Map) {
        const QPointF q = mapData(pos);
        const int b = std::max(first_, std::min(last_, (int)std::floor(q.x())));
        const int d = std::max(0, std::min(frame_.nfft - 1, (int)std::floor(q.y())));
        curB_ = b;  curD_ = d;
        const double v = cellDb(b, d);
        readout_ = QStringLiteral("R %1 km   v %2 m/s   %3")
                       .arg((b + 0.5) * frame_.dR_m / 1e3, 0, 'f', 2).arg(frame_.vel[d], 0, 'f', 1)
                       .arg(v == v ? QStringLiteral("%1 dB above noise").arg(v, 0, 'f', 1) : QStringLiteral("no cell"));
        if (const Detection *dt = detectionAt(pos)) {
            tip = cellText(binOf(dt->range_m), dt->dIdx, true);
            key = (long long)binOf(dt->range_m) * 4096 + dt->dIdx;
        } else {
            tip = cellText(b, d, true);
            key = (long long)b * 4096 + d;
        }
    } else if (a == Lower) {
        const std::vector<LowPoint> pts = lowerPoints();
        double bestD = 1e18;
        const LowPoint *best = nullptr;
        for (const LowPoint &pt : pts) {
            const double dist = std::fabs(lowScreen(pt.x, pt.db).x() - pos.x());
            if (dist < bestD) { bestD = dist; best = &pt; }
        }
        if (best && bestD < 12.0) {
            const Pin c = lowerCell(best->bin);
            tip = cellText(c.b, c.d, true);
            key = -100 - ((long long)c.b * 4096 + c.d);
        }
    }
    tip_ = tip;
    if (key != tipKey_) {
        tipKey_ = key;
        if (!tip.isEmpty()) QToolTip::showText(mapToGlobal(pos.toPoint()), tip, this);
        else QToolTip::hideText();
    }
    update();
}

void RdWidget::zoomAt(Area a, QPointF pos, double fx, double fy)
{
    if (a == Map) {
        const QPointF c = mapData(pos);
        vx0_ = c.x() - (c.x() - vx0_) / fx;  vx1_ = c.x() + (vx1_ - c.x()) / fx;
        vy0_ = c.y() - (c.y() - vy0_) / fy;  vy1_ = c.y() + (vy1_ - c.y()) / fy;
        mapAuto_ = false;
        clampMap();
        image_ = QImage();
    } else if (a == Lower) {
        const QPointF c = lowData(pos);
        lx0_ = c.x() - (c.x() - lx0_) / fx;  lx1_ = c.x() + (lx1_ - c.x()) / fx;
        ly0_ = c.y() - (c.y() - ly0_) / fy;  ly1_ = c.y() + (ly1_ - c.y()) / fy;
        lowAuto_ = false;
        clampLower();
    }
    update();
}

void RdWidget::panBy(Area a, QPointF px)
{
    if (a == Map) {
        const double dx = -px.x() / mapRect_.width() * (vx1_ - vx0_), dy = -px.y() / mapRect_.height() * (vy1_ - vy0_);
        vx0_ += dx; vx1_ += dx; vy0_ += dy; vy1_ += dy;
        mapAuto_ = false;
        clampMap();
        image_ = QImage();
    } else if (a == Lower) {
        const double dx = -px.x() / lowRect_.width() * (lx1_ - lx0_), dy = px.y() / lowRect_.height() * (ly1_ - ly0_);
        lx0_ += dx; lx1_ += dx; ly0_ += dy; ly1_ += dy;
        lowAuto_ = false;
        clampLower();
    }
    update();
}

/* click: pin the data tip of the cell, detection or lower-plot point (again: remove it);
 * a click on the map also selects the cell for the lower plot */
void RdWidget::togglePin(QPointF pos)
{
    Pin c{ -1, -1 };
    const Area a = areaAt(pos);
    if (a == Map) {
        if (const Detection *dt = detectionAt(pos)) c = Pin{ binOf(dt->range_m), dt->dIdx };
        else {
            const QPointF q = mapData(pos);
            c = Pin{ std::max(first_, std::min(last_, (int)std::floor(q.x()))),
                     std::max(0, std::min(frame_.nfft - 1, (int)std::floor(q.y()))) };
        }
        selB_ = c.b;  selD_ = c.d;
    } else if (a == Lower) {
        double bestD = 12.0;
        for (const LowPoint &pt : lowerPoints()) {
            const double dist = std::fabs(lowScreen(pt.x, pt.db).x() - pos.x());
            if (dist < bestD) { bestD = dist; c = lowerCell(pt.bin); }
        }
    }
    if (c.b < 0 || c.d < 0) return;
    auto it = std::find_if(pins_.begin(), pins_.end(), [&c](const Pin &p) { return p.b == c.b && p.d == c.d; });
    if (it != pins_.end()) pins_.erase(it);
    else pins_.push_back(c);
    update();
}

void RdWidget::mousePressEvent(QMouseEvent *e)
{
    if (e->button() != Qt::LeftButton) return;
    pressPos_ = lastPos_ = event_pos(e);
    dragArea_ = areaAt(pressPos_);
    dragging_ = false;
}

void RdWidget::mouseMoveEvent(QMouseEvent *e)
{
    const QPointF pos = event_pos(e);
    if ((e->buttons() & Qt::LeftButton) && dragArea_ != None) {
        if (!dragging_ && std::hypot(pos.x() - pressPos_.x(), pos.y() - pressPos_.y()) > 3.0) dragging_ = true;
        if (dragging_) {
            panBy(dragArea_, pos - lastPos_);
            lastPos_ = pos;
            setCursor(Qt::ClosedHandCursor);
            return;
        }
    }
    updateCursor(pos);
}

void RdWidget::mouseReleaseEvent(QMouseEvent *e)
{
    if (e->button() != Qt::LeftButton) return;
    if (!dragging_ && dragArea_ != None) togglePin(event_pos(e));
    dragging_ = false;
    dragArea_ = None;
    unsetCursor();
    updateCursor(event_pos(e));                       /* the readout and the cut at the new position */
}

void RdWidget::mouseDoubleClickEvent(QMouseEvent *e)
{
    const Area a = areaAt(event_pos(e));
    if (a == Map) mapAuto_ = true;
    else if (a == Lower) lowAuto_ = true;
    else return;
    image_ = QImage();
    fitViews();
    update();
}

void RdWidget::wheelEvent(QWheelEvent *e)
{
    const QPointF pos = event_pos(e);
    const Area a = areaAt(pos);
    if (a == None) { e->ignore(); return; }
    const double steps = e->angleDelta().y() / 120.0;
    if (steps == 0.0) return;
    const double f = std::pow(1.25, steps);
    const bool onlyX = e->modifiers() & Qt::ShiftModifier, onlyY = e->modifiers() & Qt::ControlModifier;
    zoomAt(a, pos, onlyY ? 1.0 : f, onlyX ? 1.0 : f);
    e->accept();
}

void RdWidget::keyPressEvent(QKeyEvent *e)
{
    if (!hasFrame()) { QWidget::keyPressEvent(e); return; }
    const QPointF mid = mapRect_.center(), lmid = lowRect_.center();
    const QPoint local = mapFromGlobal(QCursor::pos());
    const Area a = areaAt(local) == Lower ? Lower : Map;      /* the plot under the cursor, else the map */
    const QPointF c = a == Lower ? lmid : mid;
    const double step = 40.0;
    switch (e->key()) {
    case Qt::Key_Plus: case Qt::Key_Equal: zoomAt(a, c, 1.5, 1.5); break;
    case Qt::Key_Minus:                    zoomAt(a, c, 1.0 / 1.5, 1.0 / 1.5); break;
    case Qt::Key_Left:  panBy(a, QPointF(step, 0.0)); break;
    case Qt::Key_Right: panBy(a, QPointF(-step, 0.0)); break;
    case Qt::Key_Up:    panBy(a, QPointF(0.0, step)); break;
    case Qt::Key_Down:  panBy(a, QPointF(0.0, -step)); break;
    case Qt::Key_R: case Qt::Key_Home: resetViews(); break;
    case Qt::Key_C:     pins_.clear(); selB_ = selD_ = -1; update(); break;
    default: QWidget::keyPressEvent(e); return;
    }
}

void RdWidget::contextMenuEvent(QContextMenuEvent *e)
{
    QMenu menu(this);
    QObject::connect(menu.addAction(QStringLiteral("Reset zoom (R, double click)")), &QAction::triggered,
                     [this]() { resetViews(); });
    QObject::connect(menu.addAction(QStringLiteral("Remove data tips and selection (C)")), &QAction::triggered,
                     [this]() { pins_.clear(); selB_ = selD_ = -1; update(); });
    QAction *spec = menu.addAction(QStringLiteral("Lower plot: Doppler spectrum"));
    QAction *prof = menu.addAction(QStringLiteral("Lower plot: range profile"));
    spec->setCheckable(true);  spec->setChecked(lower_ == RdLower::Spectrum);
    prof->setCheckable(true);  prof->setChecked(lower_ == RdLower::RangeProfile);
    QObject::connect(spec, &QAction::triggered, [this]() { setLower(RdLower::Spectrum); });
    QObject::connect(prof, &QAction::triggered, [this]() { setLower(RdLower::RangeProfile); });
    menu.exec(e->globalPos());
}

void RdWidget::drawPins(QPainter &p) const
{
    const QFont font = pixel_font(this->font(), 11);
    const QFontMetrics fm(font);
    p.setFont(font);
    for (const Pin &pin : pins_) {
        const QStringList rows = cellText(pin.b, pin.d, false).split(QChar('\n'));
        int w = 0;
        for (const QString &r : rows) w = std::max(w, fm.horizontalAdvance(r));
        const double rowH = fm.height(), boxW = w + 12.0, boxH = rows.size() * rowH + 8.0;
        auto box_at = [&](QPointF s, const QRectF &area) {
            QRectF box(s.x() + 12.0, s.y() + 8.0, boxW, boxH);
            if (box.right() > width()) box.moveRight(s.x() - 12.0);
            if (box.bottom() > area.bottom() + 30.0) box.moveBottom(s.y() - 8.0);
            p.setPen(QPen(kColorPin, 1.5));
            p.setBrush(Qt::NoBrush);
            p.drawRect(QRectF(s.x() - 6.0, s.y() - 6.0, 12.0, 12.0));
            p.drawLine(s, box.topLeft());
            p.setBrush(QColor(10, 40, 10, 225));
            p.drawRect(box);
            p.setPen(Qt::white);
            for (int i = 0; i < rows.size(); ++i)
                p.drawText(QPointF(box.left() + 6.0, box.top() + 4.0 + fm.ascent() + i * rowH), rows[i]);
        };
        /* on the map, when the cell is in view */
        const QPointF s = mapScreen(pin.b + 0.5, pin.d + 0.5);
        if (mapRect_.contains(s)) box_at(s, mapRect_);
        /* on the lower plot, when the cut passes the cell */
        const bool onCut = lower_ == RdLower::Spectrum ? pin.b == spectrumBin() : pin.d == profileBin();
        const double v = cellDb(pin.b, pin.d);
        if (onCut && v == v) {
            const QPointF q = lowScreen(lowerX(lower_ == RdLower::Spectrum ? pin.d : pin.b), v);
            if (lowRect_.contains(q)) box_at(q, lowRect_);
        }
    }
}

void RdWidget::drawLower(QPainter &p, const QFontMetrics &fm)
{
    const bool spec = lower_ == RdLower::Spectrum;
    const int sb = spectrumBin(), pd = profileBin();
    p.setPen(kColorLabel);
    const QString where = curB_ >= 0 ? QStringLiteral("cursor") : (selB_ >= 0 ? QStringLiteral("selected cell")
                                                                               : QStringLiteral("strongest peak"));
    const QString title = spec
        ? QStringLiteral("Doppler spectrum at %1 km (range bin %2, %3); circles: detections")
              .arg((sb + 0.5) * frame_.dR_m / 1e3, 0, 'f', 2).arg(sb).arg(where)
        : QStringLiteral("Range profile at %1 m/s (Doppler bin %2, %3); circles: detections")
              .arg(pd >= 0 && pd < frame_.nfft ? frame_.vel[pd] : 0.0, 0, 'f', 1).arg(pd).arg(where);
    p.drawText(QPointF(lowRect_.left(), lowRect_.top() - 7), fm.elidedText(title, Qt::ElideRight, (int)lowRect_.width() + 40));
    /* grid, ticks and labels of both axes */
    const double xs = nice_step(lx1_ - lx0_, 8), ys = nice_step(ly1_ - ly0_, 5);
    for (double x = std::ceil(lx0_ / xs) * xs; x <= lx1_ + 1e-9; x += xs) {
        const double sx = lowScreen(x, ly0_).x();
        p.setPen(QPen(kColorGrid, 1.0, Qt::DotLine));
        p.drawLine(QPointF(sx, lowRect_.top()), QPointF(sx, lowRect_.bottom()));
        p.setPen(kColorLabel);
        p.drawLine(QPointF(sx, lowRect_.bottom()), QPointF(sx, lowRect_.bottom() + 4));
        const QString t = QString::number(std::fabs(x) < 1e-9 ? 0.0 : x, 'g', 5);
        p.drawText(QPointF(sx - fm.horizontalAdvance(t) / 2.0, lowRect_.bottom() + 16), t);
    }
    for (double y = std::ceil(ly0_ / ys) * ys; y <= ly1_ + 1e-9; y += ys) {
        const double sy = lowScreen(lx0_, y).y();
        p.setPen(QPen(kColorGrid, 1.0, Qt::DotLine));
        p.drawLine(QPointF(lowRect_.left(), sy), QPointF(lowRect_.right(), sy));
        p.setPen(kColorLabel);
        p.drawLine(QPointF(lowRect_.left() - 4, sy), QPointF(lowRect_.left(), sy));
        const QString t = QString::number(std::fabs(y) < 1e-9 ? 0.0 : y, 'g', 4);
        p.drawText(QPointF(lowRect_.left() - 8 - fm.horizontalAdvance(t), sy + 4), t);
    }
    const QString xl = spec ? QStringLiteral("Radial velocity [m/s]") : QStringLiteral("Range [km]");
    p.drawText(QPointF(lowRect_.center().x() - fm.horizontalAdvance(xl) / 2.0, lowRect_.bottom() + 31), xl);
    const QString yl = QStringLiteral("Power [dB above noise]");
    p.save();
    p.translate(14, lowRect_.center().y() + fm.horizontalAdvance(yl) / 2.0);
    p.rotate(-90);
    p.drawText(QPointF(0, 0), yl);
    p.restore();
    p.setPen(QPen(kColorRim, 1.0));
    p.setBrush(Qt::NoBrush);
    p.drawRect(lowRect_);
    /* the curve and the detections on it */
    const std::vector<LowPoint> pts = lowerPoints();
    p.save();
    p.setClipRect(lowRect_);
    QPainterPath path;
    for (size_t i = 0; i < pts.size(); i++) {
        const QPointF q = lowScreen(pts[i].x, pts[i].db);
        if (i == 0) path.moveTo(q); else path.lineTo(q);
    }
    p.setPen(QPen(QColor(0, 230, 255), 1.5));
    p.drawPath(path);
    p.setPen(QPen(Qt::white, 1.5));
    for (const Detection &dt : frame_.dets) {
        const int b = binOf(dt.range_m);
        if (spec ? b != sb : dt.dIdx != pd) continue;
        const double v = cellDb(b, dt.dIdx);
        if (v == v) p.drawEllipse(lowScreen(lowerX(spec ? dt.dIdx : b), v), 5.0, 5.0);
    }
    p.restore();
}

void RdWidget::paintEvent(QPaintEvent *)
{
    QPainter p(this);
    p.fillRect(rect(), Qt::black);
    layoutRects();
    const QFont small = pixel_font(font(), 11), title = pixel_font(font(), 14, true);
    const QFontMetrics fm(small);
    p.setFont(title);
    p.setPen(Qt::white);
    const QString pause = paused_ ? QStringLiteral("   [PAUSED]") : QString();
    if (!hasFrame()) {
        p.drawText(QPointF(16.0, 26.0), (frame_.valid ? QStringLiteral("Range-Doppler map: no range cell yet")
                                                      : QStringLiteral("Range-Doppler map: waiting for a CPI update of %1 ...")
                                                            .arg(beam_text(g_beam.load()))) + pause);
        return;
    }
    const int textW = std::max(50, width() - 24);              /* long lines are cut with "..." */
    p.drawText(QPointF(16.0, 22.0), QFontMetrics(title).elidedText(
                   QStringLiteral("Range-Doppler map, %1   t = %2 s   az %3   peak %4 dB above noise%5")
                       .arg(beam_text(frame_.beam)).arg(frame_.time_s, 0, 'f', 2)
                       .arg(deg_text(frame_.az_deg, 1)).arg(maxDb_, 0, 'f', 1).arg(pause), Qt::ElideRight, textW));
    p.setFont(small);
    p.setPen(kColorLabel);
    p.drawText(QPointF(16.0, 40.0), fm.elidedText(readout_.isEmpty()
                       ? QStringLiteral("wheel: zoom (Shift: range / x only, Ctrl: velocity / y only)   drag: pan   "
                                        "double click: reset   hover: data tip   click: pin, select the cut   C: remove tips")
                       : readout_, Qt::ElideRight, textW));
    const int nfft = frame_.nfft;
    const double lo = -3.0, hi = hiDb_;
    const int W = std::max(1, (int)mapRect_.width()), H = std::max(1, (int)mapRect_.height());
    if (image_.isNull() || image_.width() != W || image_.height() != H) {
        /* one pixel per pixel of the view; a pixel shows the largest of its range bins
         * (a peak of one range bin stays visible when zoomed out) */
        image_ = QImage(W, H, QImage::Format_RGB32);
        std::vector<int> b0(W), b1(W);
        const double sx = (vx1_ - vx0_) / W;
        for (int c = 0; c < W; c++) {
            const double xa = vx0_ + c * sx, xb = xa + sx;
            b0[c] = std::max(first_, std::min(last_, (int)std::floor(xa)));
            b1[c] = std::max(b0[c] + 1, std::min(last_ + 1, (int)std::ceil(xb - 1e-9)));
        }
        for (int r = 0; r < H; r++) {
            const int d = std::max(0, std::min(nfft - 1, (int)std::floor(vy0_ + (r + 0.5) / H * (vy1_ - vy0_))));
            QRgb *line = reinterpret_cast<QRgb *>(image_.scanLine(r));
            const float *row = &db_[(size_t)d * frame_.nBins];
            for (int c = 0; c < W; c++) {
                float v = NAN;
                for (int b = b0[c]; b < b1[c]; b++) if (row[b] == row[b] && !(row[b] <= v)) v = row[b];
                line[c] = v == v ? colour_db(v, lo, hi) : qRgb(0, 0, 0);
            }
        }
    }
    p.drawImage(mapRect_.topLeft(), image_);
    p.save();
    p.setClipRect(mapRect_);
    p.setRenderHint(QPainter::Antialiasing, true);
    /* MTI notch: the Doppler bins around zero that the CFAR does not test */
    const int z = nfft / 2;
    p.setPen(QPen(QColor(200, 200, 200), 1.0, Qt::DashLine));
    for (int e : { z - frame_.notchBins, z + frame_.notchBins + 1 }) {
        const double y = mapScreen(vx0_, e).y();
        p.drawLine(QPointF(mapRect_.left(), y), QPointF(mapRect_.right(), y));
    }
    /* detections of the beam */
    p.setPen(QPen(Qt::white, 1.5));
    p.setBrush(Qt::NoBrush);
    for (const Detection &dt : frame_.dets) p.drawEllipse(screenOf(dt.range_m, dt.dIdx), 6.0, 6.0);
    /* the cut of the lower plot */
    p.setPen(QPen(QColor(255, 255, 255, 120), 1.0, Qt::DotLine));
    if (lower_ == RdLower::Spectrum) {
        const double x = mapScreen(spectrumBin() + 0.5, 0).x();
        p.drawLine(QPointF(x, mapRect_.top()), QPointF(x, mapRect_.bottom()));
    } else {
        const double y = mapScreen(0, profileBin() + 0.5).y();
        p.drawLine(QPointF(mapRect_.left(), y), QPointF(mapRect_.right(), y));
    }
    p.restore();
    /* axes: range (km) below, radial velocity (m/s) left; frame */
    p.setPen(QPen(kColorRim, 1.0));
    p.drawRect(mapRect_);
    const double r0 = vx0_ * frame_.dR_m / 1e3, r1 = vx1_ * frame_.dR_m / 1e3;
    const double rs = nice_step(r1 - r0, 8);
    p.setPen(kColorLabel);
    for (double r = std::ceil(r0 / rs) * rs; r <= r1 + 1e-9; r += rs) {
        const double x = mapRect_.left() + (r - r0) / (r1 - r0) * mapRect_.width();
        p.drawLine(QPointF(x, mapRect_.bottom()), QPointF(x, mapRect_.bottom() + 4));
        const QString t = QString::number(r, 'g', 5);
        p.drawText(QPointF(x - fm.horizontalAdvance(t) / 2.0, mapRect_.bottom() + 16), t);
    }
    const QString rl = QStringLiteral("Range [km]");
    p.drawText(QPointF(mapRect_.center().x() - fm.horizontalAdvance(rl) / 2.0, mapRect_.bottom() + 31), rl);
    const double vTop = velAt(vy0_), vBot = velAt(vy1_);
    const double vs = nice_step(std::fabs(vTop - vBot), 8);
    for (double v = std::ceil(std::min(vTop, vBot) / vs) * vs; v <= std::max(vTop, vBot) + 1e-9; v += vs) {
        const double y = mapScreen(vx0_, rowOf(v)).y();
        p.drawLine(QPointF(mapRect_.left() - 4, y), QPointF(mapRect_.left(), y));
        const QString t = QString::number(std::fabs(v) < 1e-9 ? 0.0 : v, 'g', 4);
        p.drawText(QPointF(mapRect_.left() - 8 - fm.horizontalAdvance(t), y + 4), t);
    }
    const QString vl = QStringLiteral("Radial velocity [m/s]");
    p.save();
    p.translate(14, mapRect_.center().y() + fm.horizontalAdvance(vl) / 2.0);
    p.rotate(-90);
    p.drawText(QPointF(0, 0), vl);
    p.restore();
    /* colour bar */
    const QRectF bar(mapRect_.right() + 14, mapRect_.top(), 14, mapRect_.height());
    for (int i = 0; i < (int)bar.height(); i++) {
        const double x = hi - (hi - lo) * i / bar.height();
        p.setPen(QColor(colour_db(x, lo, hi)));
        p.drawLine(QPointF(bar.left(), bar.top() + i), QPointF(bar.right(), bar.top() + i));
    }
    p.setPen(kColorLabel);
    const double cs = nice_step(hi - lo, 6);
    for (double x = 0; x <= hi + 1e-9; x += cs) {
        const double y = bar.top() + (hi - x) / (hi - lo) * bar.height();
        p.drawText(QPointF(bar.right() + 4, y + 4), QString::number(x, 'f', 0));
    }
    p.drawText(QPointF(bar.right() + 4, bar.bottom() + 14), QStringLiteral("dB"));
    /* the lower plot: Doppler spectrum or range profile through the cut */
    p.setRenderHint(QPainter::Antialiasing, true);
    drawLower(p, fm);
    drawPins(p);
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
        mode_->addItem(QStringLiteral("Video (DDC, elements avg.)"));
        mode_->addItem(QStringLiteral("DBF out (selected beam)"));
        mode_->addItem(QStringLiteral("MF out (selected beam)"));
        v->addWidget(mode_);
        QObject::connect(mode_, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
                         [this](int i) { ppi_->setDisplayMode((DisplayMode)i); ppi_->setFocus(); });
        QLabel *beamLabel = new QLabel(QStringLiteral("Beam after DBF (DBF out, MF out, range-Doppler)"), panel);
        beamLabel->setWordWrap(true);
        v->addWidget(beamLabel);
        beam_ = new QComboBox(panel);
        for (int b = 0; b < (int)beam_elevations().size(); b++) {
            const double e = beam_elevations()[b];
            beam_->addItem(std::isfinite(e) ? QStringLiteral("%1   (el %2)").arg(b + 1).arg(deg_text(e, 1))
                                            : QStringLiteral("%1   (no DBF)").arg(b + 1));
        }
        v->addWidget(beam_);
        QObject::connect(beam_, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), [this](int i) {
            ppi_->setBeam(i);
            rd_->setFrame(RdFrame());                    /* the map of the new beam comes with the next update */
            ppi_->setFocus();
        });
        tracks_ = new QCheckBox(QStringLiteral("Tracks (T)"), panel);
        tracks_->setChecked(ppi_->showTracks());
        v->addWidget(tracks_);
        QObject::connect(tracks_, &QCheckBox::toggled, [this](bool on) { ppi_->setShowTracks(on); ppi_->setFocus(); });
        QLabel *lowerLabel = new QLabel(QStringLiteral("Range-Doppler tab: plot below the map"), panel);
        lowerLabel->setWordWrap(true);
        v->addWidget(lowerLabel);
        lower_ = new QComboBox(panel);
        lower_->addItem(QStringLiteral("Doppler spectrum (range bin)"));
        lower_->addItem(QStringLiteral("Range profile (Doppler bin)"));
        v->addWidget(lower_);
        QObject::connect(lower_, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
                         [this](int i) { rd_->setLower(i == 1 ? RdLower::RangeProfile : RdLower::Spectrum); });
        pause_ = new QPushButton(QStringLiteral("Pause display (Space)"), panel);
        pause_->setCheckable(true);
        v->addWidget(pause_);
        QObject::connect(pause_, &QPushButton::toggled, [this](bool on) { setPaused(on); });

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
        timing_->setMinimumHeight(QFontMetrics(mono).lineSpacing() * 15 + 12);  /* 15 rows, stable layout */
        timing_->setAlignment(Qt::AlignLeft | Qt::AlignTop);
        v->addWidget(timing_);
        QLabel *help = new QLabel(QStringLiteral("Wheel: zoom   Drag: pan\nDouble click: reset\nHover: data tip   "
                                                 "Click: pin\nKeys: + - arrows R C D T S Space\nRight click: more options\n"
                                                 "(PPI and range-Doppler plots)"),
                                  panel);
        help->setStyleSheet(QStringLiteral("color: #909090;"));
        v->addWidget(help);

        /* tabs: the PPI; the range-Doppler map of the selected beam with a small PPI.
         * The one PPI widget moves between the tabs; the panel serves both. */
        tabs_ = new QTabWidget(this);
        tabs_->setStyleSheet(QStringLiteral("QTabWidget::pane { border: 0; } QTabBar::tab { background: #202020; "
                                            "color: #d0d0d0; padding: 5px 16px; } QTabBar::tab:selected { "
                                            "background: #383838; color: white; }"));
        QWidget *page1 = new QWidget(), *page2 = new QWidget();
        /* the pages only (by name): an unscoped "background: black" would also reach the
         * data tips (QToolTip of the PPI), drawn black on black */
        page1->setObjectName(QStringLiteral("ppiPage"));
        page2->setObjectName(QStringLiteral("rdPage"));
        page1->setStyleSheet(QStringLiteral("QWidget#ppiPage { background: black; }"));
        page2->setStyleSheet(QStringLiteral("QWidget#rdPage { background: black; }"));
        pageLayout1_ = new QVBoxLayout(page1);
        pageLayout1_->setContentsMargins(0, 0, 0, 0);
        pageLayout2_ = new QHBoxLayout(page2);
        pageLayout2_->setContentsMargins(0, 0, 0, 0);
        pageLayout2_->setSpacing(4);
        rd_ = new RdWidget();
        pageLayout2_->addWidget(rd_, 3);
        pageLayout1_->addWidget(ppi_, 1);
        tabs_->addTab(page1, QStringLiteral("PPI"));
        tabs_->addTab(page2, QStringLiteral("Range-Doppler (selected beam)"));
        QObject::connect(tabs_, &QTabWidget::currentChanged, [this](int i) { showTab(i); });

        QHBoxLayout *top = new QHBoxLayout(this);
        top->setContentsMargins(0, 0, 0, 0);
        top->setSpacing(0);
        top->addWidget(tabs_, 1);
        top->addWidget(panel);
        ppi_->setFocus();
    }

    PpiWidget &ppi() { return *ppi_; }
    RdWidget &rd() { return *rd_; }
    QComboBox *lowerBox() { return lower_; }
protected:
    void keyPressEvent(QKeyEvent *e) override       /* Space in the range-Doppler tab (the PPI has its own) */
    {
        if (e->key() == Qt::Key_Space) setPaused(!ppi_->paused());
        else QWidget::keyPressEvent(e);
    }
public:
    QPushButton *pauseButton() { return pause_; }
    /* pause / resume the display: the PPI and the range-Doppler map take no new data
     * (the processing goes on; the plots that came meanwhile are shown on resume) */
    void setPaused(bool on)
    {
        ppi_->setPaused(on);
        rd_->setPaused(on);
        const bool b = pause_->blockSignals(true);
        pause_->setChecked(on);
        pause_->setText(on ? QStringLiteral("Resume display (Space)") : QStringLiteral("Pause display (Space)"));
        pause_->setStyleSheet(on ? QStringLiteral("background: #5a4a10; color: white;") : QString());
        pause_->blockSignals(b);
    }
    QTabWidget *tabs() { return tabs_; }
    QComboBox *beamBox() { return beam_; }
    /* the PPI in the first tab (whole page), or small next to the range-Doppler map */
    void showTab(int i)
    {
        if (i == 1) pageLayout2_->addWidget(ppi_, 2);
        else pageLayout1_->addWidget(ppi_, 1);
        ppi_->setMinimumSize(i == 1 ? 200 : 300, i == 1 ? 200 : 300);
        g_rdWanted.store(i == 1);
        if (tabs_->currentIndex() != i) tabs_->setCurrentIndex(i);
        ppi_->setFocus();
    }
    /* screenshot: the PPI in the first tab, the whole page (map and PPI) in the second */
    bool saveSnapshot(const QString &path)
    {
        if (tabs_->currentIndex() == 0) return ppi_->saveSnapshot(path);
        const QPixmap pm = tabs_->currentWidget()->grab();
        if (pm.save(path, "PNG")) return true;
        std::fprintf(stderr, "PPI: cannot save the snapshot to %s\n", path.toLocal8Bit().constData());
        return false;
    }
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
        if (pause_->isChecked() != ppi_->paused()) setPaused(ppi_->paused());   /* Space or the context menu */
        if (lower_->currentIndex() != (int)rd_->lower()) {         /* changed in the map's context menu */
            const bool b = lower_->blockSignals(true);
            lower_->setCurrentIndex((int)rd_->lower());
            lower_->blockSignals(b);
        }
        if (tracks_->isChecked() != ppi_->showTracks()) {          /* changed with T or the context menu */
            const bool b = tracks_->blockSignals(true);
            tracks_->setChecked(ppi_->showTracks());
            tracks_->blockSignals(b);
        }
        showTiming();
        {                                                   /* newest range-Doppler map (second tab) */
            RdFrame f;
            bool fresh = false;
            {
                Shared &s = shared();
                std::lock_guard<std::mutex> lock(s.mutex);
                if (s.rdSerial != rdSerial_ && tabs_->currentIndex() == 1 && !ppi_->paused()) {   /* paused: kept */
                    rdSerial_ = s.rdSerial;
                    if (s.rd.beam == ppi_->beam()) { f = s.rd; fresh = true; }
                }
            }
            if (fresh) rd_->setFrame(f);
        }
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
        /* the mean over the PRIs since the last shown values (at least 100),
         * not since the start: a change (e.g. a lower CFAR threshold) shows
         * at once instead of creeping into a long-term mean */
        if (t.pris < timingPrev_.pris) timingPrev_ = ProcessorTiming{}, guiPrev_ = 0.0;   /* new run */
        const uint64_t nPri = t.pris - timingPrev_.pris;
        if (nPri < 100 && timingPrev_.pris > 0) return;
        timingSerial_ = serial;
        const bool newName = detName_ != QString::fromStdString(det);
        detName_ = QString::fromStdString(det);
        if (newName) cfarSerial_ = 0;                          /* rewrite the CFAR line with the name */
        const ProcessorTiming &o = timingPrev_;
        const double k = 1e3 / (double)nPri;                   /* ms in the interval -> us per PRI */
        /* processing thread: the blocks, the console output and the data handed
         * to the PPI (all within the PRI: Total); the drawing runs in the GUI
         * thread, in parallel with the processing, so it is not part of Total.
         * Labels of at most 15 characters, so the value column stays visible. */
        const std::pair<QString, double> rows[] = {
            { QStringLiteral("DBF"), t.dbf_ms - o.dbf_ms }, { QStringLiteral("Matched filter"), t.mf_ms - o.mf_ms },
            { QStringLiteral("MTI (clutter)"), t.mti_ms - o.mti_ms },
            { QStringLiteral("Doppler FFT+NCI"), t.doppler_ms - o.doppler_ms },
            { detName_, t.cfar_ms - o.cfar_ms }, { QStringLiteral("Peak interp."), t.interp_ms - o.interp_ms },
            { QStringLiteral("Plot extractor"), t.extract_ms - o.extract_ms },
            { QString::fromStdString(trk), t.track_ms - o.track_ms },
            { QStringLiteral("Console output"), t.console_ms - o.console_ms },
            { QStringLiteral("Display (post)"), t.display_ms - o.display_ms } };
        QStringList lines;
        double sum = 0.0;
        for (const auto &r : rows) {
            lines << timing_row(r.first, r.second * k);
            sum += r.second * k;
        }
        lines << timing_row(QStringLiteral("Total"), sum);
        lines << timing_row(QStringLiteral("PPI draw (GUI)"), (ppi_->guiMs() - guiPrev_) * k);
        lines << QStringLiteral(" (GUI thread, in parallel:");
        lines << QStringLiteral("  not part of Total)");
        lines << QStringLiteral("last %1 of %2 PRIs").arg((qulonglong)nPri).arg((qulonglong)t.pris);
        timingPrev_ = t;
        guiPrev_ = ppi_->guiMs();
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
    QComboBox             *mode_, *beam_, *lower_;
    QCheckBox             *tracks_;
    QPushButton           *pause_;
    QTabWidget            *tabs_;
    QVBoxLayout           *pageLayout1_;
    QHBoxLayout           *pageLayout2_;
    RdWidget              *rd_;
    uint64_t               rdSerial_ = 0;
    QLabel                *timing_;
    uint64_t               timingSerial_ = 0;
    ProcessorTiming        timingPrev_ = {};    /* accumulated times at the start of the shown interval */
    double                 guiPrev_ = 0.0;
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

extern "C" int ppi_video_wanted(void)
{
    return g_videoWanted.load();
}

extern "C" void ppi_post_video(const float *row, int nBins, int nAvg, double dR_m, double az_deg, double time_s)
{
    if (!row || nBins <= 0 || !(dR_m > 0.0) || !std::isfinite(az_deg)) return;
    VideoRow r;
    r.nAvg = nAvg > 0 ? nAvg : 1;
    r.az_deg = az_deg;
    r.time_s = time_s;
    r.dR_m = dR_m;
    r.v.assign(row, row + nBins);
    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (s.newVideo.size() < kMaxPendingVideo) s.newVideo.push_back(std::move(r));
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
static std::atomic<int> g_startTab(0);
static std::atomic<int> g_startRdPlot(0);

extern "C" void ppi_set_beams(int nBeams, const double *beamEl_deg)
{
    std::vector<double> &el = beam_elevations();
    el.assign(nBeams > 0 ? (size_t)nBeams : 1, NAN);
    for (int b = 0; b < nBeams && beamEl_deg; b++) el[b] = beamEl_deg[b];
}

extern "C" void ppi_set_beam(int beam)
{
    g_beam.store(beam >= 0 && beam < (int)beam_elevations().size() ? beam : 0);
}

extern "C" void ppi_set_tab(int tab)
{
    g_startTab.store(tab == 1);
}

extern "C" void ppi_set_rd_plot(int plot)
{
    g_startRdPlot.store(plot == 1);
}

extern "C" int ppi_profile_wanted(int *beam)
{
    if (beam) *beam = g_beam.load();
    return g_profileKind.load();
}

extern "C" void ppi_post_profile(int kind, const float *row, int nBins, double dR_m, double az_deg, double time_s)
{
    if (!row || nBins <= 0 || !(dR_m > 0.0) || !std::isfinite(az_deg)) return;
    VideoRow r;
    r.kind = kind;
    r.beam = g_beam.load();
    r.az_deg = az_deg;
    r.time_s = time_s;
    r.dR_m = dR_m;
    r.v.assign(row, row + nBins);
    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (s.newVideo.size() < kMaxPendingVideo) s.newVideo.push_back(std::move(r));
}

extern "C" int ppi_rdmap_wanted(int *beam)
{
    if (beam) *beam = g_beam.load();
    return g_rdWanted.load();
}

extern "C" void ppi_post_rdmap(const RdMap *map)
{
    if (!map || map->nBins <= 0 || map->nfft <= 0 || !map->power || !map->vel_mps) return;
    RdFrame f;
    f.valid = true;
    f.beam = map->beam;
    f.nBins = map->nBins;
    f.nfft = map->nfft;
    f.notchBins = map->notchBins;
    f.dR_m = map->dR_m;
    f.noise = map->noise;
    f.time_s = map->time_s;
    f.az_deg = map->az_deg;
    f.power.assign(map->power, map->power + (size_t)map->nBins * map->nfft);
    f.vel.assign(map->vel_mps, map->vel_mps + map->nfft);
    if (map->dets && map->nDets > 0) f.dets.assign(map->dets, map->dets + map->nDets);
    Shared &s = shared();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.rd = std::move(f);
    s.rdSerial++;
}

extern "C" void ppi_set_display_mode(int mode)
{
    if (mode >= PPI_SHOW_PLOTS && mode <= PPI_SHOW_MF) g_startMode.store(mode);
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
    widget.setDisplayMode((DisplayMode)g_startMode.load());    /* (also when the index did not change) */
    window.beamBox()->setCurrentIndex(g_beam.load());
    window.showTab(g_startTab.load());
    window.lowerBox()->setCurrentIndex(g_startRdPlot.load());
    window.showMaximized();                                    /* the whole screen at the start */

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
            QCoreApplication::exit(window.saveSnapshot(QString::fromLocal8Bit(snapshotPath)) ? 0 : 1);
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

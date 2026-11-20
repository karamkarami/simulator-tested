/*
 * test_ppi_interactive.cpp - interaction with the Qt PPI.
 *
 * Includes ppi_qt.cpp to reach the widget, sends real mouse, wheel and key
 * events (off-screen) and checks: hover data tips, cursor readout, wheel zoom
 * around the cursor, drag pan, click to pin a data tip, double click reset,
 * keys, show/hide and SNR filters, plot persistence options, pause; the CFAR
 * detections before the plot extractor (display modes, arcs, data tips, key
 * D) and the window panel (display mode box, CFAR threshold Apply/Default
 * through ppi_take_cfar_request / ppi_post_cfar_threshold, screenshot,
 * mean block times, CFAR method name); tracks (triangle drawn around the
 * plot, data tip, Tracks checkbox and key T, deletion).
 * Saves snapshots (argv[1]: folder, optional) of the zoomed PPI with a pinned
 * tip and of the detections.
 */
#include "../src/ppi_qt.cpp"

#include <QElapsedTimer>
#include <QEventLoop>
#include <QThread>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

namespace {

int failures = 0;

void check(bool ok, const char *what)
{
    std::printf("  %-60s %s\n", what, ok ? "ok" : "FAIL");
    if (!ok) failures++;
}

void send_mouse(QWidget &w, QEvent::Type type, QPointF pos, Qt::MouseButton button)
{
    const Qt::MouseButtons buttons = (type == QEvent::MouseButtonRelease || type == QEvent::MouseMove)
                                         ? (type == QEvent::MouseMove && button != Qt::NoButton ? Qt::MouseButtons(button) : Qt::MouseButtons(Qt::NoButton))
                                         : Qt::MouseButtons(button);
    QMouseEvent e(type, pos, QPointF(w.mapToGlobal(pos.toPoint())), button, buttons, Qt::NoModifier);
    QCoreApplication::sendEvent(&w, &e);
}

void send_wheel(QWidget &w, QPointF pos, int angle, Qt::KeyboardModifiers mods = Qt::NoModifier)
{
    QWheelEvent e(pos, w.mapToGlobal(pos.toPoint()), QPoint(), QPoint(0, angle), Qt::NoButton, mods,
                  Qt::NoScrollPhase, false);
    QCoreApplication::sendEvent(&w, &e);
}

void send_key(QWidget &w, int key, const char *text)
{
    QKeyEvent e(QEvent::KeyPress, key, Qt::NoModifier, QString::fromLatin1(text));
    QCoreApplication::sendEvent(&w, &e);
}

Detection make_plot(double range_m, double az_deg, double el_deg, double vr, double snr, int pulse, uint32_t number)
{
    Detection d;
    std::memset(&d, 0, sizeof d);
    d.range_m = range_m;  d.az_deg = az_deg;  d.el_deg = el_deg;  d.vr_mps = vr;  d.snr_dB = snr;
    d.pulse = pulse;  d.hits = 70;  d.beam = 4;  d.rIdx = 123;  d.dIdx = 45;  d.time_s = 1.5;
    d.pulseNumber = number;
    return d;
}

QPointF ground(const Detection &d)
{
    const double g = d.range_m * (std::isfinite(d.el_deg) ? std::cos(d.el_deg * kDeg) : 1.0);
    return QPointF(g * std::sin(d.az_deg * kDeg), g * std::cos(d.az_deg * kDeg));
}

} /* namespace */

int main(int argc, char **argv)
{
    if (qEnvironmentVariableIsEmpty("QT_QPA_PLATFORM")) qputenv("QT_QPA_PLATFORM", "offscreen");
    const std::string folder = argc > 1 ? argv[1] : ".";
    QApplication app(argc, argv);
    reset_shared();
    PpiWidget w(60e3, nullptr);
    w.resize(800, 800);
    w.show();

    const Detection a = make_plot(10e3, 30.0, NAN, 12.5, 25.0, PULSE_SHORT, 101);
    const Detection b = make_plot(40e3, 200.0, 3.0, -7.25, 14.0, PULSE_LONG, 102);
    const Detection c = make_plot(30e3, 95.0, 5.0, 3.0, 30.0, PULSE_LONG, 103);
    ppi_post_antenna(10.0, 0.1, 1);
    ppi_post_plot(&a);
    ppi_post_antenna(100.0, 2.6, 1);
    ppi_post_plot(&c);
    ppi_post_antenna(210.0, 5.6, 1);
    ppi_post_plot(&b);
    ppi_post_antenna(215.0, 5.8, 1);
    w.refresh();
    check(w.plotCount() == 3, "three plots on the display");

    /* data tips */
    const QPointF ga = ground(a), gc = ground(c);
    QPointF sa = w.screenOf(ga.x(), ga.y());
    QString tip = w.tipText(sa + QPointF(3, 2));
    check(tip.contains("short pulse") && tip.contains("10.000 km") && tip.contains("12.50 m/s") &&
          tip.contains("no DBF") && tip.contains("Pulse number 101"), "hover data tip shows the plot values");
    check(w.tipText(sa + QPointF(30, 30)).isEmpty(), "no data tip away from the plots");
    send_mouse(w, QEvent::MouseMove, sa, Qt::NoButton);
    check(w.cursorText().contains("km") && w.cursorText().contains(QChar(0x00B0)), "cursor readout (range, azimuth)");

    /* wheel zoom around the cursor */
    const QPointF sc = w.screenOf(gc.x(), gc.y());
    send_wheel(w, sc, 120);
    const QPointF sc2 = w.screenOf(gc.x(), gc.y());
    check(std::fabs(w.zoom() - 1.25) < 1e-9 && std::hypot(sc2.x() - sc.x(), sc2.y() - sc.y()) < 0.5,
          "wheel zooms x1.25, the point under the cursor stays");
    for (int i = 0; i < 8; ++i) send_wheel(w, sc, 120);
    check(w.tipText(w.screenOf(gc.x(), gc.y())).contains("30.000 km"), "data tip works in the zoomed view");

    /* drag = pan */
    const QPointF v0 = w.viewCentre_m();
    const double k = w.screenOf(1000.0, 0.0).x() - w.screenOf(0.0, 0.0).x();   /* pixels per km */
    send_mouse(w, QEvent::MouseButtonPress, QPointF(400, 400), Qt::LeftButton);
    send_mouse(w, QEvent::MouseMove, QPointF(460, 400), Qt::LeftButton);
    send_mouse(w, QEvent::MouseButtonRelease, QPointF(460, 400), Qt::LeftButton);
    const QPointF v1 = w.viewCentre_m();
    check(std::fabs((v0.x() - v1.x()) - 60.0 / k * 1000.0) < 1.0 && std::fabs(v0.y() - v1.y()) < 1.0 &&
          w.pinnedCount() == 0, "drag pans the view (no pin)");

    /* click on a plot pins its data tip, a second click removes it */
    const QPointF sc3 = w.screenOf(gc.x(), gc.y());
    send_mouse(w, QEvent::MouseButtonPress, sc3, Qt::LeftButton);
    send_mouse(w, QEvent::MouseButtonRelease, sc3, Qt::LeftButton);
    check(w.pinnedCount() == 1, "click pins the data tip");
    if (!folder.empty()) w.saveSnapshot(QString::fromStdString(folder + "/ppi_zoom_pinned.png"));
    send_mouse(w, QEvent::MouseButtonPress, sc3, Qt::LeftButton);
    send_mouse(w, QEvent::MouseButtonRelease, sc3, Qt::LeftButton);
    check(w.pinnedCount() == 0, "second click removes it");

    /* double click and keys */
    send_mouse(w, QEvent::MouseButtonDblClick, QPointF(400, 400), Qt::LeftButton);
    check(w.zoom() == 1.0 && w.viewCentre_m() == QPointF(0, 0), "double click resets the view");
    send_key(w, Qt::Key_Plus, "+");
    check(std::fabs(w.zoom() - 1.5) < 1e-9, "key + zooms in");
    send_key(w, Qt::Key_Left, "");
    check(w.viewCentre_m().x() < 0.0, "arrow key pans");
    send_key(w, Qt::Key_R, "r");
    check(w.zoom() == 1.0, "key R resets the view");
    send_wheel(w, QPointF(400, 400), -120);
    check(w.zoom() == 1.0, "no zoom out beyond the full view");

    /* filters */
    w.setShowShort(false);
    sa = w.screenOf(ga.x(), ga.y());
    check(w.tipText(sa).isEmpty() && w.visibleCount() == 2, "short-pulse plots hidden (no data tip)");
    w.setShowShort(true);
    w.setMinSnr(20.0);
    check(w.visibleCount() == 2, "minimum SNR 20 dB hides the 14 dB plot");
    w.setMinSnr(-1e9);

    /* persistence: the antenna passes plot a (30 deg) again in scan 2 */
    w.setHistory(History::All);
    ppi_post_antenna(359.0, 9.9, 1);
    ppi_post_antenna(40.0, 11.2, 2);
    w.refresh();
    check(w.plotCount() == 3, "'all' keeps the plots");
    w.setHistory(History::UntilSwept);
    check(w.plotCount() == 2, "default: a plot disappears when the antenna passes it again");

    /* pause */
    w.setPaused(true);
    const Detection d = make_plot(20e3, 45.0, 2.0, 1.0, 20.0, PULSE_SHORT, 104);
    ppi_post_antenna(50.0, 11.3, 2);
    ppi_post_plot(&d);
    w.refresh();
    check(w.plotCount() == 2, "paused display does not take new plots");
    w.setPaused(false);
    w.refresh();
    check(w.plotCount() == 3, "resumed display shows them");
    w.clearPlots();
    check(w.plotCount() == 0, "clear plots");

    if (!folder.empty()) w.saveSnapshot(QString::fromStdString(folder + "/ppi_cleared.png"));

    /* CFAR detections before the plot extractor: a target seen in 20 CPI updates */
    w.setHistory(History::All);
    ppi_post_antenna(60.0, 12.0, 2);
    const double arc = 6.0 * 1e-3 * 8.0 * 6.0;           /* 6 rpm = 36 deg/s, 1 ms PRI, every 8 pulses */
    for (int i = 0; i < 20; i++) {
        Detection dt = make_plot(25e3, 70.0 + i * arc, 2.0, -3.0, 18.0, PULSE_LONG, 200 + i);
        dt.hits = 0;
        ppi_post_detections(&dt, 1, arc);
    }
    const Detection dsh = make_plot(12e3, 75.0, NAN, 5.0, 11.0, PULSE_SHORT, 300);
    ppi_post_detections(&dsh, 1, arc);
    ppi_post_plot(&c);
    w.refresh();
    check(w.detectionCount() == 21 && w.plotCount() == 1, "detections and plots are stored separately");
    check(w.displayMode() == DisplayMode::Plots && w.visibleDetectionCount() == 0 && w.visibleCount() == 1,
          "default: plots only");
    const Detection mid = make_plot(25e3, 70.0 + 10 * arc, 2.0, -3.0, 18.0, PULSE_LONG, 210);
    const QPointF gm = ground(mid);
    check(w.tipText(w.screenOf(gm.x(), gm.y())).isEmpty(), "no detection data tip in the plots mode");
    w.setDisplayMode(DisplayMode::Detections);
    check(w.visibleDetectionCount() == 21 && w.visibleCount() == 0, "detections mode: detections only");
    for (int i = 0; i < 6; ++i) send_wheel(w, w.screenOf(gm.x(), gm.y()), 120);
    tip = w.tipText(w.screenOf(gm.x(), gm.y()));
    check(tip.contains("CFAR detection, long pulse") && tip.contains("25.000 km") && !tip.contains("CPI hits"),
          "hover data tip of a detection");
    if (!folder.empty()) w.saveSnapshot(QString::fromStdString(folder + "/ppi_detections.png"));
    {   /* the arcs are drawn: pixels of the detection colour along the arc */
        const QImage img = w.grab().toImage();
        const QPointF sm = w.screenOf(gm.x(), gm.y());
        const QColor px = img.pixelColor(sm.toPoint());
        check(px.red() > 100 && px.green() > 70 && px.blue() < 60, "the arc is drawn in the detection colour");
    }
    w.setMinSnr(15.0);
    check(w.visibleDetectionCount() == 20, "SNR filter hides the 11 dB detection");
    w.setMinSnr(-1e9);
    send_key(w, Qt::Key_D, "d");
    check(w.displayMode() == DisplayMode::Both && w.visibleDetectionCount() == 21 && w.visibleCount() == 1,
          "key D: plots and detections");
    send_key(w, Qt::Key_D, "d");
    check(w.displayMode() == DisplayMode::Video && w.visibleDetectionCount() == 0 && w.visibleCount() == 0,
          "key D again: the raw video (no plots, no detections)");
    send_key(w, Qt::Key_D, "d");
    check(w.displayMode() == DisplayMode::DbfOut && w.visibleCount() == 0, "key D again: DBF out of the selected beam");
    send_key(w, Qt::Key_D, "d");
    check(w.displayMode() == DisplayMode::MfOut && w.visibleCount() == 0, "key D again: MF out of the selected beam");
    send_key(w, Qt::Key_D, "d");
    check(w.displayMode() == DisplayMode::Plots, "key D again: plots");
    w.resetView();

    /* window: display mode box, CFAR threshold, screenshot */
    const double beamEl[3] = { 1.5, 4.0, 7.5 };     /* the beam selector: 3 beams */
    ppi_set_beams(3, beamEl);
    PpiWindow win(60e3, nullptr);
    win.resize(1030, 800);
    win.show();
    win.refresh();
    check(!win.applyButton()->isEnabled() && win.cfarInfo().contains("waiting"), "CFAR field waits for the processing");
    ppi_post_cfar_threshold(12.79, 1.0e-8, 12.79, 9.0);
    win.refresh();
    check(win.applyButton()->isEnabled() && std::fabs(win.cfarBox()->value() - 12.79) < 1e-9 &&
          win.cfarInfo().contains(": 12.79 dB") && win.cfarInfo().contains("1.0e-08"),
          "the threshold reported by the processing is shown");
    double req = 0.0;
    check(ppi_take_cfar_request(&req) == 0, "no request before Apply");
    win.cfarBox()->setValue(17.5);
    win.applyButton()->click();
    check(ppi_take_cfar_request(&req) == 1 && req == 17.5 && win.cfarInfo().contains("requested"),
          "Apply sends the new threshold to the processing");
    check(ppi_take_cfar_request(&req) == 0, "the request is taken once");
    win.cfarBox()->setValue(5.0);
    check(std::fabs(win.cfarBox()->value() - 9.0) < 1e-9, "the threshold field stops at the threshold of cfg.cfar.maxPfa");
    ppi_post_cfar_threshold(17.5, 2e-20, 12.79, 9.0);
    win.refresh();
    check(win.cfarInfo().contains(": 17.50 dB"), "the applied threshold is shown");
    win.defaultButton()->click();
    check(ppi_take_cfar_request(&req) == 1 && std::fabs(req - 12.79) < 1e-12, "Default requests the design threshold");
    win.modeBox()->setCurrentIndex(2);
    check(win.ppi().displayMode() == DisplayMode::Both, "display mode box switches the PPI");
    send_key(win.ppi(), Qt::Key_D, "d");
    win.refresh();
    check(win.modeBox()->currentIndex() == 3 && ppi_video_wanted() == 1, "the box follows a mode change by key D (video)");
    for (int k = 4; k <= 5; k++) {
        send_key(win.ppi(), Qt::Key_D, "d");
        win.refresh();
        int bm = -1;
        check(win.modeBox()->currentIndex() == k && ppi_video_wanted() == 0 &&
              ppi_profile_wanted(&bm) == (k == 4 ? PROFILE_DBF : PROFILE_MF) && bm == 0,
              k == 4 ? "key D: DBF out, the profile of beam 1 wanted" : "key D: MF out, the profile of beam 1 wanted");
    }
    send_key(win.ppi(), Qt::Key_D, "d");
    win.refresh();
    check(win.modeBox()->currentIndex() == 0 && ppi_video_wanted() == 0 && ppi_profile_wanted(nullptr) == PROFILE_NONE,
          "key D from MF out back to the plots");
    {   /* the raw outputs offered (ppi_set_raw_outputs): SLC out only with the SLC, DBF out only with the DBF */
        ppi_set_raw_outputs(1, 1);
        PpiWindow ws(60e3, nullptr);
        check(ws.modeBox()->count() == 7 && ws.modeBox()->itemText(5).startsWith("SLC out") &&
              ws.modeBox()->itemText(6).startsWith("MF out"), "with the SLC: the box lists SLC out between DBF out and MF out");
        ws.modeBox()->setCurrentIndex(5);
        int bm = -1;
        check(ws.ppi().displayMode() == DisplayMode::SlcOut && ppi_profile_wanted(&bm) == PROFILE_SLC,
              "the box selects SLC out: the profile of the SLC output wanted");
        send_key(ws.ppi(), Qt::Key_D, "d");
        ws.refresh();
        check(ws.ppi().displayMode() == DisplayMode::MfOut && ws.modeBox()->currentIndex() == 6, "key D: from SLC out to MF out");
        ws.modeBox()->setCurrentIndex(4);
        send_key(ws.ppi(), Qt::Key_D, "d");
        check(ws.ppi().displayMode() == DisplayMode::SlcOut, "key D: from DBF out to SLC out");
        ws.modeBox()->setCurrentIndex(0);
        ppi_set_raw_outputs(0, 1);
        PpiWindow wn(60e3, nullptr);
        check(wn.modeBox()->count() == 6 && wn.modeBox()->findData((int)DisplayMode::DbfOut) < 0 &&
              wn.modeBox()->findData((int)DisplayMode::SlcOut) == 4, "without the DBF (one main element): no DBF out, SLC out");
        wn.modeBox()->setCurrentIndex(0);
        ppi_set_raw_outputs(1, 0);
        PpiWindow wd(60e3, nullptr);
        check(wd.modeBox()->count() == 6 && wd.modeBox()->findData((int)DisplayMode::SlcOut) < 0,
              "without the SLC: no SLC out (the default)");
        wd.modeBox()->setCurrentIndex(0);
    }
    {   /* the data tip in the window (in its tab page) is readable: text unlike its background */
        const Detection t = make_plot(25e3, 45.0, 2.0, 50.0, 20.0, PULSE_LONG, 777);
        ppi_post_antenna(45.0, 7.0, 1);
        ppi_post_plot(&t);
        ppi_post_antenna(50.0, 7.1, 1);
        win.refresh();
        const QPointF gt = ground(t);
        send_mouse(win.ppi(), QEvent::MouseMove, win.ppi().screenOf(gt.x(), gt.y()) + QPointF(2, 1), Qt::NoButton);
        QApplication::processEvents();
        QWidget *tipW = nullptr;
        for (QWidget *tw : QApplication::allWidgets())
            if (tw->inherits("QTipLabel")) tipW = tw;
        int lo = 255, hi = 0;
        if (tipW) {
            const QImage im = tipW->grab().toImage();
            for (int y = 0; y < im.height(); y++)
                for (int x = 0; x < im.width(); x++) {
                    const int g = qGray(im.pixel(x, y));
                    lo = std::min(lo, g);
                    hi = std::max(hi, g);
                }
            tipW->grab().save(QString::fromStdString(folder + "/ppi_window_datatip.png"));
        }
        check(tipW && QToolTip::text().contains("Pulse number 777") && hi - lo > 100,
              "data tip in the window: shown, with its text readable (not black on black)");
        QToolTip::hideText();
    }
    /* block times in the panel, CFAR method name */
    ProcessorTiming tm;
    std::memset(&tm, 0, sizeof tm);
    tm.pris = 1000;
    tm.dbf_ms = 120.0;  tm.mf_ms = 340.0;  tm.mti_ms = 56.0;  tm.doppler_ms = 980.0;
    tm.cfar_ms = 760.0; tm.extract_ms = 15.0; tm.track_ms = 4.0;
    tm.console_ms = 30.0; tm.display_ms = 50.0;
    tm.interp_ms = 7.0;
    ppi_post_timing(&tm, "CFAR 1-D range", "Tracking 3-D CV with a long name");
    win.refresh();
    const QString timing = win.timingText();
    auto row = [](const char *name, double us) { return QStringLiteral("%1%2").arg(QString(name), -15).arg(us, 9, 'f', 1); };
    check(timing.contains("Doppler FFT+NCI") && timing.contains("980.0") && timing.contains(row("CFAR 1D range", 760.0)) &&
          timing.contains(row("Peak interp.", 7.0)) && timing.contains(row("Tracking 3D CV.", 4.0)) &&
          timing.contains(row("Total", 2362.0)), "mean block times per PRI in the panel (labels shortened)");
    check(timing.contains(row("Console output", 30.0)) && timing.contains(row("Display (post)", 50.0)) &&
          timing.contains("PPI draw (GUI)"), "console, display and GUI drawing times in the panel");
    check(timing_label("CFAR OS 1-D Doppler") == "CFAR OS 1D Dop" && timing_label("CFAR GO 1-D range") == "CFAR GO 1D rng" &&
          timing_label("CFAR CA 2-D") == "CFAR CA 2D" && timing_label("Track 3-D PMBM") == "Track 3D PMBM",
          "the labels of every CFAR detector, window and tracker fit in 15 characters");
    {
        int longest = 0;
        for (const QString &l : timing.split(QChar('\n')).mid(0, 12)) longest = std::max(longest, (int)l.length());
        check(longest <= 24, "every timing row fits in 24 characters (the value is not pushed out)");
    }
    win.grab();                                          /* draws the PPI in the GUI thread */
    tm.pris = 1200;  tm.cfar_ms += 400.0;                /* 200 more PRIs, the CFAR block slower */
    ppi_post_timing(&tm, "CFAR 1-D range", "Track 3D CV");
    win.refresh();
    const QString timing2 = win.timingText();
    const int iGui = timing2.indexOf("PPI draw (GUI)");
    check(iGui >= 0 && timing2.mid(iGui + 15, 9).trimmed().toDouble() > 0.0, "the drawing time of the GUI is counted");
    check(timing2.contains(row("CFAR 1D range", 2000.0)) && timing2.contains(row("Total", 2000.0)) &&
          timing2.contains("last 200 of 1200 PRIs") && timing2.contains("not in Total") &&
          timing2.contains("Logging (copy)") && timing2.contains("Log writer"),
          "the times are the means over the PRIs since the last update (GUI drawing not in Total)");
    check(win.cfarInfo().contains("CFAR 1-D range: 17.50 dB"), "the CFAR line names the method");

    /* tracks: hollow triangles around their plot, checkbox, deletion */
    PpiWidget &pw = win.ppi();
    pw.setHistory(History::All);
    ppi_post_plot(&c);                                   /* plot c at 30 km, 95 deg */
    Track tk;
    std::memset(&tk, 0, sizeof tk);
    tk.id = 7;  tk.confirmed = 1;  tk.dim = 3;  tk.time_s = 12.0;  tk.hits = 5;
    tk.pos_m[0] = gc.x() + 40.0;  tk.pos_m[1] = gc.y();  tk.pos_m[2] = 2600.0;  tk.alt_m = 2606.0;
    tk.vel_mps[0] = -100.0;  tk.vel_mps[1] = 50.0;
    ppi_post_track(&tk);
    Track tt = tk;
    tt.id = 8;  tt.confirmed = 0;  tt.hits = 1;  tt.pos_m[0] = 5e3;  tt.pos_m[1] = 5e3;
    ppi_post_track(&tt);
    win.refresh();
    check(pw.trackCount() == 2 && pw.confirmedCount() == 1, "a confirmed and a tentative track");
    for (int i = 0; i < 6; ++i) send_wheel(pw, pw.screenOf(gc.x(), gc.y()), 120);
    const QPointF st = pw.screenOf(gc.x(), gc.y());
    tip = pw.tipText(st);
    check(tip.contains("Track 7 (confirmed, 3-D EKF)") && tip.contains("Height: 2606 m") && tip.contains("Speed: 111.8 m/s"),
          "hover over a track (wins over its plot)");
    if (!folder.empty()) pw.saveSnapshot(QString::fromStdString(folder + "/ppi_track.png"));
    {   /* the triangle is drawn around the plot: magenta pixels above the plot */
        const QImage img = pw.grab().toImage();
        bool magenta = false;
        for (int dy = 6; dy <= 14 && !magenta; ++dy) {
            const QColor px = img.pixelColor(QPoint((int)std::lround(st.x()), (int)std::lround(st.y()) - dy));
            magenta = px.red() > 150 && px.blue() > 150 && px.green() < 140;
        }
        check(magenta, "the track triangle is drawn around the plot");
    }
    {   /* a PMBM track: existence and the IMM model probabilities in its data tip */
        Track tp = tk;
        tp.method = TRACK_METHOD_PMBM;
        tp.existence = 0.973;
        tp.nModes = 3;
        tp.modes[0] = PMBM_CV;  tp.modes[1] = PMBM_CA;  tp.modes[2] = PMBM_CT;
        tp.modeProb[0] = 0.05;  tp.modeProb[1] = 0.15;  tp.modeProb[2] = 0.80;
        ppi_post_track(&tp);
        win.refresh();
        const QString tipP = pw.tipText(st);
        check(tipP.contains("Track 7 (confirmed, 3-D PMBM)") && tipP.contains("Existence: 0.973") &&
              tipP.contains("Models: CV 0.05 CA 0.15 CT 0.80"), "data tip of a PMBM track: existence, model probabilities");
        ppi_post_track(&tk);                             /* back to the EKF track */
        win.refresh();
    }
    win.tracksBox()->setChecked(false);
    check(!pw.showTracks() && !pw.tipText(st).contains("Track"), "Tracks checkbox hides the tracks");
    send_key(pw, Qt::Key_T, "t");
    win.refresh();
    check(pw.showTracks() && win.tracksBox()->isChecked(), "key T shows them again (checkbox follows)");
    tk.deleted = 1;
    ppi_post_track(&tk);
    win.refresh();
    check(pw.trackCount() == 1 && pw.confirmedCount() == 0, "a deleted track disappears");
    pw.resetView();

    /* raw video: rows of the worker (video.c) for 30 ... 70 deg, an echo at 30 km for 45 ... 55 deg */
    {
        win.modeBox()->setCurrentIndex(3);
        check(pw.displayMode() == DisplayMode::Video && ppi_video_wanted() == 1, "the box selects the video");
        const double dR = 74.948;
        const int nBins = (int)std::ceil(60e3 / dR);
        std::vector<float> row(nBins);
        uint32_t seed = 1;
        auto noise = [&seed]() {                              /* one element: exponential power of mean 800 */
            seed = seed * 1664525u + 1013904223u;
            return (float)(-800.0 * std::log(((seed >> 8) + 0.5) / 16777216.0));
        };
        for (double az = 30.0; az < 70.0; az += 0.036) {
            for (int i = 0; i < nBins; i++) row[i] = i < 55 ? 0.0f : noise();   /* no data inside 4.1 km */
            if (az > 45.0 && az < 55.0)
                for (int i = 400; i < 404; i++) row[i] = 40000.0f;
            ppi_post_video(row.data(), nBins, 1, dR, az, 20.0 + az / 36.0);
        }
        ppi_post_antenna(70.0, 22.0, 2);
        win.refresh();
        check(std::fabs(p10_offset_dB(1) - 9.7732) < 1e-3 && std::fabs(p10_offset_dB(48) - 0.8616) < 1e-3,
              "noise floor offset of the 10th percentile: one sample 9.77 dB, mean of 48 0.86 dB (gammaincinv)");
        check(pw.videoBinsFilled() == 160 && std::fabs(pw.videoNoise_dB() - 10.0 * std::log10(800.0)) < 0.5,
              "video of 30 ... 70 deg (160 azimuth bins), noise floor of one element = its mean power");
        const double r = 30.1e3;
        check(pw.videoIndexAt(r * std::sin(50 * kDeg), r * std::cos(50 * kDeg)) > 80 &&
              pw.videoIndexAt(20e3 * std::sin(40 * kDeg), 20e3 * std::cos(40 * kDeg)) < 40 &&   /* < 5 dB above noise */
              pw.videoIndexAt(20e3 * std::sin(150 * kDeg), 20e3 * std::cos(150 * kDeg)) == 0,
              "echo bright, noise dark, no video where the antenna did not pass");
        const QImage img = pw.grab().toImage();
        const QPointF se = pw.screenOf(r * std::sin(50 * kDeg), r * std::cos(50 * kDeg));
        int best = 0;
        for (int dy = -2; dy <= 2; dy++)
            for (int dx = -2; dx <= 2; dx++) {
                const QColor px = img.pixelColor(QPoint((int)se.x() + dx, (int)se.y() + dy));
                best = std::max(best, px.red() + px.green() + px.blue());
            }
        check(best > 250 && pw.visibleCount() == 0, "the echo is drawn on the PPI; plots are hidden in video mode");
        if (!folder.empty()) pw.saveSnapshot(QString::fromStdString(folder + "/ppi_video.png"));
        pw.clearPlots();
        check(pw.videoBinsFilled() == 0, "Clear display also clears the video");
        win.modeBox()->setCurrentIndex(0);
        check(ppi_video_wanted() == 0, "back to the plots: no video wanted");
    }

    /* video noise floor over a mountain sector: 48 elements (mean of 48 exponential powers), the
     * rows of 180 ... 225 deg have clutter 30 dB above the noise in 95% of their bins (near
     * clutter, long-pulse smear and the mountain range); the noise floor and so the colour of the
     * noise must not move while the antenna crosses it */
    {
        win.modeBox()->setCurrentIndex(3);
        const double dR = 74.948;
        const int nBins = (int)std::ceil(80e3 / dR);
        std::vector<float> row(nBins);
        uint32_t seed = 7;
        std::vector<float> pool(65536);                       /* means of 48 exponential powers of mean 500 */
        for (float &x : pool) {
            double sum = 0.0;
            for (int k = 0; k < 48; k++) {
                seed = seed * 1664525u + 1013904223u;
                sum += -500.0 * std::log(((seed >> 8) + 0.5) / 16777216.0);
            }
            x = (float)(sum / 48.0);
        }
        auto noise48 = [&seed, &pool]() {
            seed = seed * 1664525u + 1013904223u;
            return pool[seed >> 16];
        };
        double lo = 1e9, hi = -1e9, before = NAN;
        int n = 0;
        for (double az = 60.0; az < 300.0; az += 0.036) {   /* 1 ms PRI at 6 rpm */
            const bool mountain = az > 180.0 && az < 225.0;
            for (int i = 0; i < nBins; i++)
                row[i] = i < 55 ? 0.0f : noise48() * (mountain && i < 55 + (int)(0.95 * (nBins - 55)) ? 1000.0f : 1.0f);
            ppi_post_video(row.data(), nBins, 48, dR, az, 100.0 + az / 36.0);
            if (++n % 50) continue;
            win.refresh();
            const double v = pw.videoNoise_dB();
            if (az > 150.0 && az < 170.0) before = v;
            if (az > 150.0) { lo = std::min(lo, v); hi = std::max(hi, v); }
        }
        char what[200];
        std::snprintf(what, sizeof what, "video noise floor over a mountain sector (95%% clutter in 45 deg): %.1f dB, "
                      "%.2f dB of variation", before, hi - lo);
        check(std::fabs(before - 10.0 * std::log10(500.0)) < 0.5 && hi - lo < 0.5, what);
        pw.clearPlots();
        win.modeBox()->setCurrentIndex(0);
    }

    /* DBF out / MF out of the selected beam */
    {
        check(win.beamBox()->count() == 3 && win.beamBox()->itemText(1).contains("4.0"), "beam selector: 3 beams with elevations");
        win.beamBox()->setCurrentIndex(2);
        int bm = -1;
        check(pw.beam() == 2 && ppi_profile_wanted(&bm) == PROFILE_NONE && bm == 2, "beam 3 selected, no profile in plot mode");
        win.modeBox()->setCurrentIndex(4);
        check(ppi_profile_wanted(&bm) == PROFILE_DBF && bm == 2 && ppi_video_wanted() == 0, "DBF out: the profile of beam 3 wanted");
        const double dR = 74.948;
        const int nBins = (int)std::ceil(60e3 / dR);
        std::vector<float> row(nBins);
        uint32_t seed = 3;
        for (double az = 100.0; az < 120.0; az += 0.036) {
            for (int i = 0; i < nBins; i++) {                     /* exponential noise of mean 500 (one sample) */
                seed = seed * 1664525u + 1013904223u;
                row[i] = i < 55 ? 0.0f : (float)(-500.0 * std::log(((seed >> 8) + 0.5) / 16777216.0));
            }
            row[300] = 50000.0f;                                  /* an echo at 22.5 km */
            ppi_post_profile(PROFILE_DBF, row.data(), nBins, dR, az, 30.0 + az / 36.0);
            ppi_post_profile(PROFILE_MF, row.data(), nBins, dR, az, 30.0 + az / 36.0);   /* other kind: dropped */
            ppi_post_video(row.data(), nBins, 1, dR, az + 180.0, 30.0 + az / 36.0);           /* video: dropped */
        }
        win.refresh();
        const double r = (300 + 0.5) * dR;
        check(pw.videoBinsFilled() == 80 && pw.videoIndexAt(r * std::sin(110 * kDeg), r * std::cos(110 * kDeg)) > 100 &&
              pw.videoIndexAt(20e3 * std::sin(290 * kDeg), 20e3 * std::cos(290 * kDeg)) == 0,
              "DBF out: only the DBF rows (100 ... 120 deg), the echo bright, no video drawn");
        check(std::fabs(pw.videoNoise_dB() - 10.0 * std::log10(500.0)) < 1.0,
              "DBF out: noise floor from one sample per bin (10th percentile + 9.77 dB)");
        const QImage img = pw.grab().toImage();
        check(!img.isNull(), "the PPI draws the DBF output");
        if (!folder.empty()) pw.saveSnapshot(QString::fromStdString(folder + "/ppi_dbf_out.png"));
        win.modeBox()->setCurrentIndex(5);
        check(pw.videoBinsFilled() == 0 && ppi_profile_wanted(&bm) == PROFILE_MF, "MF out: the DBF data cleared, MF wanted");
        ppi_post_profile(PROFILE_MF, row.data(), nBins, dR, 50.0, 31.0);
        win.beamBox()->setCurrentIndex(1);                        /* posted for beam 3: no longer shown */
        win.refresh();
        check(pw.videoBinsFilled() == 0 && pw.beam() == 1, "a profile of the previous beam is not shown");
        ppi_post_profile(PROFILE_MF, row.data(), nBins, dR, 50.0, 31.1);
        win.refresh();
        check(pw.videoBinsFilled() == 1, "MF out of beam 2 shown");
        win.modeBox()->setCurrentIndex(0);
        check(ppi_profile_wanted(nullptr) == PROFILE_NONE, "back to the plots: no profile wanted");
    }

    /* the range-Doppler tab: the map of the selected beam next to a small PPI */
    {
        int bm = -1;
        check(win.tabs()->count() == 3 && win.tabs()->currentIndex() == 0 && ppi_rdmap_wanted(&bm) == 1 &&
              pw.parentWidget() == win.tabs()->widget(0), "three tabs; the map is computed also in the first");
        win.resize(1400, 800);
        win.tabs()->setCurrentIndex(1);
        QCoreApplication::processEvents();
        check(ppi_rdmap_wanted(&bm) == 1 && bm == 1 && pw.parentWidget() == win.tabs()->widget(1) &&
              !win.rd().hasFrame(), "second tab: map of beam 2 wanted, the PPI moved next to it");
        {   /* before the first map: the empty map and plot below it with grid and labels; the selector
             * of the lower plot on top of that plot */
            RdWidget &r0 = win.rd();
            const QImage e = r0.grab().toImage();
            int grid = 0;                                     /* grey grid / rim pixels inside the map */
            const QRectF mr = r0.mapRect();
            for (int x = (int)mr.left() + 2; x < (int)mr.right() - 1; x++)
                for (int y = (int)mr.top() + 2; y < (int)mr.bottom() - 1; y++) {
                    const QColor q = e.pixelColor(x, y);
                    grid += q.red() > 20 && std::abs(q.red() - q.blue()) < 20;   /* (fill: 6, 6, 6) */
                }
            const QRect lb = win.lowerBox()->geometry();
            check(grid > 20 && win.lowerBox()->parentWidget() == &r0 && win.lowerBox()->isVisible() &&
                  lb.bottom() < r0.lowerRect().top() && lb.top() > r0.mapRect().bottom(),
                  "empty range-Doppler tab: grid and axes drawn, the plot selector on top of the lower plot");
            if (!folder.empty()) win.grab().save(QString::fromStdString(folder + "/ppi_rd_empty.png"));
        }
        const int nB = 800, nfft = 64, pb = 400, pd = 20;
        const double dR = 74.948;
        std::vector<float> P((size_t)nB * nfft);
        std::vector<double> vel(nfft);
        uint32_t seed = 7;
        for (size_t i = 0; i < P.size(); i++) {
            seed = seed * 1664525u + 1013904223u;
            P[i] = (i % nB) < 100 ? 0.0f : 0.5f + (float)(seed >> 8) / 16777216.0f;   /* no cells below bin 100 */
        }
        for (int d = pd - 2; d <= pd + 2; d++) P[(size_t)d * nB + pb] = (float)(1000.0 * std::exp(-0.5 * (d - pd) * (d - pd)));
        for (int d = 0; d < nfft; d++) vel[d] = -(d - nfft / 2) * 6.9;
        Detection det;
        std::memset(&det, 0, sizeof det);
        det.range_m = (pb + 0.5) * dR;  det.dIdx = pd;  det.beam = 1;  det.power = 1000.0;
        RdMap m;
        std::memset(&m, 0, sizeof m);
        m.beam = 1;  m.nBins = nB;  m.nfft = nfft;  m.dR_m = dR;  m.power = P.data();  m.vel_mps = vel.data();
        m.notchBins = 2;  m.noise = 1.0;  m.time_s = 40.0;  m.az_deg = 123.0;  m.dets = &det;  m.nDets = 1;
        ppi_post_rdmap(&m);
        win.refresh();
        RdWidget &rdw = win.rd();
        check(rdw.hasFrame() && std::fabs(rdw.maxDb() - 30.0) < 0.01 && rdw.spectrumBin() == pb,
              "map taken (peak 30 dB above noise), the spectrum at the strongest detection");
        const QImage gi = win.grab().toImage();
        const QPointF pk = rdw.mapToGlobal(rdw.screenOf(det.range_m, pd).toPoint()) - win.mapToGlobal(QPoint(0, 0));
        QColor c(0, 0, 0);                                        /* (a peak of one range bin is one pixel wide) */
        for (int dy = -2; dy <= 2; dy++)
            for (int dx = -2; dx <= 2; dx++) {
                const QColor q = gi.pixelColor(pk.toPoint() + QPoint(dx, dy));
                if (q.red() + q.green() + q.blue() > c.red() + c.green() + c.blue()) c = q;
            }
        const QPointF bg = rdw.mapToGlobal(rdw.screenOf(50e3, 50).toPoint()) - win.mapToGlobal(QPoint(0, 0));
        const QColor cb = gi.pixelColor(bg.toPoint());
        check(c.red() + c.green() + c.blue() > 400 && cb.red() + cb.green() + cb.blue() < 200,
              "the peak is drawn bright, the noise dark");
        rdw.hoverAt(rdw.screenOf(det.range_m, pd));
        check(rdw.readout().contains("30.0 dB above noise") && rdw.spectrumBin() == pb, "hover: the cell's values, its spectrum");
        /* data tips, zoom, pan, pins on the map and on the plot below it */
        check(rdw.tipText().contains("CFAR detection") && rdw.tipText().contains("30.0 dB above noise") &&
              QToolTip::text() == rdw.tipText(), "map: data tip of the detection under the cursor");
        rdw.hoverAt(rdw.screenOf(30e3, 50));
        check(rdw.tipText().contains("Range-Doppler cell") && rdw.tipText().contains("Doppler bin 50"),
              "map: data tip of a cell");
        const QRectF full = rdw.mapView();
        const QPointF at = rdw.screenOf(det.range_m, pd);
        send_wheel(rdw, at, 3 * 120);
        const QPointF at2 = rdw.screenOf(det.range_m, pd);
        check(rdw.mapView().width() < 0.6 * full.width() && rdw.mapView().height() < 0.6 * full.height() &&
              std::hypot(at2.x() - at.x(), at2.y() - at.y()) < 1.5, "map: wheel zooms both axes around the cursor");
        send_wheel(rdw, at, 120, Qt::ShiftModifier);
        const QRectF z1 = rdw.mapView();
        check(z1.height() > 0.99 * rdw.mapView().height() && z1.width() < full.width(), "map: Shift + wheel zooms the range only");
        const QRectF before = rdw.mapView();
        send_mouse(rdw, QEvent::MouseButtonPress, at, Qt::LeftButton);
        send_mouse(rdw, QEvent::MouseMove, at + QPointF(60, 0), Qt::LeftButton);
        send_mouse(rdw, QEvent::MouseMove, at + QPointF(120, 30), Qt::LeftButton);
        send_mouse(rdw, QEvent::MouseButtonRelease, at + QPointF(120, 30), Qt::LeftButton);
        const QPointF at3 = rdw.screenOf(det.range_m, pd);
        check(rdw.mapView().left() < before.left() && rdw.pinCount() == 0 &&
              std::fabs(at3.x() - at.x() - 120) < 1.5 && std::fabs(at3.y() - at.y() - 30) < 1.5,
              "map: drag pans (the peak follows the mouse), no pin");
        send_mouse(rdw, QEvent::MouseButtonPress, at3, Qt::LeftButton);
        send_mouse(rdw, QEvent::MouseButtonRelease, at3, Qt::LeftButton);
        check(rdw.pinCount() == 1, "map: click pins the data tip");
        if (!folder.empty()) rdw.grab().save(QString::fromStdString(folder + "/ppi_rd_zoom_pin.png"));
        RdFrame same = rdw.frame();
        same.time_s = 40.5;
        rdw.setFrame(same);
        check(rdw.pinCount() == 1 && std::fabs(rdw.mapView().left() - before.left()) > 1e-9 &&
              rdw.mapView().width() < full.width(), "a new map of the beam keeps the zoom and the pinned tip");
        send_mouse(rdw, QEvent::MouseButtonDblClick, at3, Qt::LeftButton);
        check(rdw.mapView() == full, "map: double click resets the zoom");
        send_mouse(rdw, QEvent::MouseButtonPress, rdw.screenOf(det.range_m, pd), Qt::LeftButton);
        send_mouse(rdw, QEvent::MouseButtonRelease, rdw.screenOf(det.range_m, pd), Qt::LeftButton);
        check(rdw.pinCount() == 0, "map: a second click removes the pinned tip");
        rdw.hoverAt(QPointF(-5, -5));                             /* the cursor leaves the map */
        check(rdw.spectrumBin() == pb, "the selected cell (the click) stays the cut of the lower plot");
        /* the plot below the map: axis labels, zoom, data tip; the range profile */
        const QPointF lp = rdw.lowerScreenOf(pd, 30.0);
        check(rdw.lowerRect().contains(lp), "spectrum: the peak inside the plot (room above it)");
        rdw.hoverAt(lp);
        check(rdw.tipText().contains("CFAR detection") && rdw.spectrumBin() == pb, "spectrum: data tip of the point under the cursor");
        const QRectF lfull = rdw.lowerView();
        send_wheel(rdw, lp, 2 * 120);
        check(rdw.lowerView().width() < 0.7 * lfull.width() && rdw.lowerView().height() < 0.7 * lfull.height(),
              "spectrum: wheel zooms around the cursor");
        send_mouse(rdw, QEvent::MouseButtonDblClick, rdw.lowerRect().center(), Qt::LeftButton);
        check(rdw.lowerView() == lfull, "spectrum: double click resets the zoom");
        win.lowerBox()->setCurrentIndex(1);
        rdw.hoverAt(QPointF(-5, -5));
        const QPointF rp = rdw.lowerScreenOf(pb, 30.0);
        rdw.hoverAt(rp);
        check(rdw.lower() == RdLower::RangeProfile && rdw.profileBin() == pd && rdw.lowerRect().contains(rp) &&
              rdw.tipText().contains("Range: 30.0") && rdw.tipText().contains("CFAR detection"),
              "panel box: the range profile through the detection's Doppler bin, its data tip");
        if (!folder.empty()) win.grab().save(QString::fromStdString(folder + "/ppi_rd_range_profile.png"));
        rdw.hoverAt(rdw.screenOf(20e3, 10));
        check(rdw.profileBin() == 10, "range profile: follows the cursor's Doppler bin");
        win.lowerBox()->setCurrentIndex(0);
        rdw.hoverAt(QPointF(-5, -5));
        /* pause: no new map and no new plots until resumed */
        win.pauseButton()->click();
        check(pw.paused() && win.pauseButton()->text().startsWith("Resume"), "pause button: the display paused");
        m.time_s = 42.0;
        ppi_post_rdmap(&m);
        const Detection late = make_plot(35e3, 300.0, 3.0, 10.0, 20.0, PULSE_LONG, 900);
        const size_t nPlots = pw.plotCount();
        ppi_post_plot(&late);
        win.refresh();
        check(rdw.frame().time_s == 40.5 && pw.plotCount() == nPlots, "paused: no new map, no new plot");
        win.pauseButton()->click();
        win.refresh();
        check(!pw.paused() && rdw.frame().time_s == 42.0 && pw.plotCount() == nPlots + 1 &&
              win.pauseButton()->text().startsWith("Pause"), "resume: the newest map and the plots of the pause shown");
        send_key(pw, Qt::Key_Space, " ");
        win.refresh();
        check(pw.paused() && win.pauseButton()->isChecked(), "the button follows Space in the PPI");
        send_key(pw, Qt::Key_Space, " ");
        win.refresh();
        m.beam = 2;  m.time_s = 41.0;
        ppi_post_rdmap(&m);                                       /* a map of another beam is not shown */
        win.refresh();
        check(rdw.frame().time_s == 42.0, "a map of another beam is not taken");
        if (!folder.empty()) {
            const QString shot = QString::fromStdString(folder + "/ppi_rd_tab.png");
            check(win.saveSnapshot(shot) && QImage(shot).width() == win.width() && QImage(shot).height() == win.height(),
                  "screenshot of the second tab: the whole window (map, plot below it, small PPI, panel)");
        }
        win.tabs()->setCurrentIndex(0);
        check(ppi_rdmap_wanted(nullptr) == 1 && pw.parentWidget() == win.tabs()->widget(0),
              "back to the PPI tab (the map is still computed)");
    }

    /* the tab Controls: logging (folder, name, blocks, start / stop), without the panel */
    {
        win.tabs()->setCurrentIndex(2);
        QCoreApplication::processEvents();
        check(!win.panel()->isVisible() && win.logButton()->text() == "Start logging",
              "Controls tab: the panel of the other tabs is hidden");
        win.logFolderEdit()->setText("/tmp/logs_test");
        win.logNameEdit()->setText("run1");
        for (int i = 0; i < 5; i++) win.logBlockBox(i)->setChecked(i != 3);
        win.logButton()->click();
        int start = -1, mask = 0;
        char folder2[256], name[64];
        const int got = ppi_take_log_request(&start, folder2, sizeof folder2, name, sizeof name, &mask);
        check(got == 1 && start == 1 && std::string(folder2) == "/tmp/logs_test" && std::string(name) == "run1" &&
              mask == (LOG_DDC | LOG_DBF | LOG_MF | LOG_DETECTION) && win.logButton()->text() == "Stop logging",
              "Start logging: the request with the folder, name and blocks");
        LoggerStats st;
        std::memset(&st, 0, sizeof st);
        st.pris = 12;  st.files = 48;  st.bytes = 3000000;  st.seconds = 65.0;
        ppi_post_log_status(1, "/tmp/logs_test/run1", LOG_DDC, &st, nullptr);
        win.refresh();
        check(pw.logging() && win.logStatus().contains("logging to /tmp/logs_test/run1") && win.logStatus().contains("48 files") &&
              !win.logBlockBox(0)->isEnabled(), "the state of the log; the blocks are locked while it runs");
        win.tabs()->setCurrentIndex(0);
        QCoreApplication::processEvents();
        {
            const QImage img = pw.grab().toImage();      /* the red REC mark in the top right corner */
            int red = 0;
            for (int x = pw.width() - 160; x < pw.width(); x++)
                for (int y = 4; y < 36; y++) {
                    const QColor q = img.pixelColor(x, y);
                    red += q.red() > 100 && q.green() < 60 && q.blue() < 60;
                }
            check(win.panel()->isVisible() && red > 20, "PPI: the REC mark while logging, the panel back in the PPI tab");
            if (!folder.empty()) win.grab().save(QString::fromStdString(folder + "/ppi_logging.png"));
        }
        win.tabs()->setCurrentIndex(2);
        win.logButton()->click();                         /* stop */
        check(ppi_take_log_request(&start, folder2, sizeof folder2, name, sizeof name, &mask) == 1 && start == 0,
              "Stop logging: the request");
        ppi_post_log_status(0, "/tmp/logs_test/run1", 0, &st, nullptr);
        win.refresh();
        check(!pw.logging() && win.logButton()->text() == "Start logging" && win.logStatus().contains("log ended") &&
              win.logBlockBox(0)->isEnabled(), "log ended: the button and the blocks back, no REC mark");
        if (!folder.empty()) win.grab().save(QString::fromStdString(folder + "/ppi_controls.png"));

        /* screen record: the window (with the mouse pointer) and the sound to <log>/screen_record.mp4 (ffmpeg;
         * the sound here from ffmpeg's test source instead of a microphone) */
        const QString ffmpeg = QStandardPaths::findExecutable("ffmpeg");
        if (ffmpeg.isEmpty()) {
            std::printf("  (no ffmpeg on PATH: the screen recording is not tested)\n");
        } else {
            ppi_set_screen_record(nullptr, "lavfi:sine=frequency=440", 10);
            const QString dir = QDir::tempPath() + "/ppi_screen_test";
            QDir(dir).removeRecursively();
            QDir().mkpath(dir);
            for (int i = 0; i < 5; i++) win.logBlockBox(i)->setChecked(false);
            win.screenRecordBox()->setChecked(true);
            win.logButton()->click();
            check(ppi_take_log_request(&start, folder2, sizeof folder2, name, sizeof name, &mask) == 1 && start == 1 &&
                  mask == LOG_SCREEN, "Start logging with only the screen record");
            ppi_post_log_status(1, dir.toUtf8().constData(), LOG_SCREEN, &st, nullptr);
            win.refresh();
            const bool rec = win.recording();
            QElapsedTimer el;
            el.start();
            while (el.elapsed() < 2500) { QCoreApplication::processEvents(QEventLoop::AllEvents, 20); QThread::msleep(5); }
            win.refresh();
            const bool shown = win.logStatus().contains("screen: recording");
            win.logButton()->click();                     /* stop: the video ends at once */
            ppi_take_log_request(&start, folder2, sizeof folder2, name, sizeof name, &mask);
            ppi_post_log_status(0, dir.toUtf8().constData(), 0, &st, nullptr);
            win.refresh();
            QProcess probe;
            probe.start(ffmpeg.left(ffmpeg.size() - 6) + "ffprobe",
                        { "-v", "error", "-show_entries", "stream=codec_type,width,height,duration", "-of", "csv=p=0",
                          dir + "/screen_record.mp4" });
            probe.waitForFinished(10000);
            const QString info = QString::fromLatin1(probe.readAllStandardOutput()).trimmed();
            const QStringList lines = info.split('\n');
            double vDur = 0, aDur = 0;
            for (const QString &l : lines) {
                const QStringList f = l.split(',');
                if (f.size() >= 4 && f[0] == "video") vDur = f[3].toDouble();
                if (f.size() >= 2 && f[0] == "audio") aDur = f.last().toDouble();
            }
            check(rec && shown && !win.recording() && win.logStatus().contains("screen recording saved") &&
                  vDur > 1.2 && vDur < 4.0 && std::fabs(aDur - vDur) < 0.3,   /* (ffmpeg starts up to ~1 s late on a busy VM) */
                  QStringLiteral("screen_record.mp4: video %1 s, sound %2 s (%3)").arg(vDur).arg(aDur)
                      .arg(info.simplified()).toStdString().c_str());
            win.screenRecordBox()->setChecked(false);
        }
        win.tabs()->setCurrentIndex(0);
    }

    if (!folder.empty()) {
        const QString shot = QString::fromStdString(folder + "/ppi_window_screenshot.png");
        check(win.ppi().saveSnapshot(shot) && QImage(shot).width() == win.width() && QImage(shot).height() == win.height(),
              "screenshot (button, key S): the whole window with the panel");
        const QIcon icon = ppi_icon();
        check(!icon.isNull() && icon.pixmap(64, 64).width() == 64, "the program icon is in the resources");
        check(QGuiApplication::desktopFileName().isEmpty() || QGuiApplication::desktopFileName() == "radar_receiver",
              "application id");
#if !defined(_WIN32) && !defined(__APPLE__)
        {                                                /* desktop entry and icon file (Linux desktops) */
            const QString home = QString::fromStdString(folder + "/desktop_home");
            QDir(home).removeRecursively();
            QString written;
            const bool okEntry = ppi_install_desktop_entry(home, &written);
            QFile e(home + "/applications/radar_receiver.desktop");
            const bool hasText = e.open(QIODevice::ReadOnly) && e.readAll().contains("Icon=" + home.toUtf8() +
                                                                                    "/icons/hicolor/256x256/apps/radar_receiver.png");
            const QImage img(home + "/icons/hicolor/256x256/apps/radar_receiver.png");
            QString again;
            ppi_install_desktop_entry(home, &again);
            check(okEntry && hasText && img.width() == 256 && !written.isEmpty() && again.isEmpty(),
                  "desktop entry radar_receiver.desktop with the icon file (written once)");
        }
#endif
        QCoreApplication::processEvents();               /* layout after the text changes */
        win.grab().save(QString::fromStdString(folder + "/ppi_window.png"));
    }
    /* Radar Settings of the tab Controls with a simulated board (board_control.c, board_link.c) */
    {
        char err[1024];
        uint8_t regs[256] = { 0 };
        regs[0x01] = 1; regs[0x02] = 20; regs[0x03] = 20; regs[0x28] = 1;
        for (int k = 0x29; k <= 0x30; k++) regs[k] = 1;
        BoardSim *sim = board_sim_open(50231, "127.0.0.1", 50232, regs, 0, err, sizeof err);
        BoardInfo bi;
        const char *sets[3] = { "cfg.board.host = '127.0.0.1';", "cfg.board.port = 50231;", "cfg.board.replyPort = 50232;" };
        const int okInfo = board_info_load(RADAR_SOURCE_DIR "/../radar_config.m", sets, 3, &bi, err, sizeof err) == 0;
        check(sim && okInfo, "Radar Settings: the simulated board and the modes of radar_config.m");
        if (sim && okInfo) {
            ModeSettings cur = board_mode_info(&bi, 0)->defaults;
            BoardControl *bc = board_control_start(&bi, 0, &cur, err, sizeof err);
            ppi_set_board_control(bc);
            PpiWindow bw(60e3, nullptr);
            bw.resize(1300, 1000);
            bw.show();
            bw.tabs()->setCurrentIndex(2);
            auto settle = [&]() {                         /* the board answers until the control is idle */
                BoardState st;
                for (int i = 0; i < 400; i++) {
                    board_sim_poll(sim, 5, nullptr);
                    board_control_state(bc, &st);
                    if (!st.busy) break;
                }
                bw.refresh();
                QCoreApplication::processEvents();
            };
            board_control_request(bc, BOARD_REQ_CHECK, 0, 0, nullptr);
            settle();
            QComboBox *code = bw.boardCombo(BOARD_CODE), *stagger = bw.boardCombo(BOARD_STAGGER);
            check(code && code->count() == 4 && code->isEnabled() && !stagger->isEnabled() &&
                  bw.boardCombo(BOARD_MODE)->count() == 3 && bw.boardCombo(BOARD_STC_J2)->count() == 4 &&
                  bw.boardCarrier(0)->value() == 170 && bw.boardStatus().contains("holds the settings"),
                  "Radar Settings: the lists of mode 0, the check of the board at the start");
            code->setCurrentIndex(2);
            emit code->activated(2);                      /* the user chooses phase code 3 */
            settle();
            check(board_sim_register(sim, 0x01) == 3 && bw.boardStatus().contains("Phase Code = 3 set"),
                  "Radar Settings: a phase code is written, read back and shown");
            RadarConfig next;
            int mode = -1;
            check(board_control_take_config(bc, &next, &mode, nullptr) == 1 && next.codeIndex[0] == 3,
                  "Radar Settings: the receiver configuration of the new phase code");
            radar_config_free(&next);
            QComboBox *modeBox = bw.boardCombo(BOARD_MODE);
            modeBox->setCurrentIndex(1);
            emit modeBox->activated(1);                   /* mode 1 (LFM: no phase code) */
            settle();
            check(board_sim_register(sim, 0x00) == 1 && !bw.boardCombo(BOARD_CODE)->isEnabled(),
                  "Radar Settings: mode 1, its phase code dropdown is off");
            while (board_control_take_config(bc, &next, &mode, nullptr)) radar_config_free(&next);
            bw.boardCarrier(1)->setValue(176);
            bw.boardCarrierApply(1)->click();
            settle();
            check(board_sim_register(sim, 0x03) == 26, "Radar Settings: F2 = 176 MHz is sent as 26");
            while (board_control_take_config(bc, &next, &mode, nullptr)) radar_config_free(&next);
            if (!folder.empty()) bw.grab().save(QString::fromStdString(folder + "/ppi_radar_settings.png"));
            ppi_set_board_control(nullptr);
            board_control_stop(bc);
        }
        board_sim_close(sim);
    }
    std::printf(failures ? "%d FAILURE(S)\n" : "PASS\n", failures);
    return failures ? 1 : 0;
}

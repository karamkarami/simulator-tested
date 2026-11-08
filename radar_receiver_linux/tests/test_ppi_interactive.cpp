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

void send_wheel(QWidget &w, QPointF pos, int angle)
{
    QWheelEvent e(pos, w.mapToGlobal(pos.toPoint()), QPoint(), QPoint(0, angle), Qt::NoButton, Qt::NoModifier,
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
    ppi_post_cfar_threshold(12.79, 1.0e-8, 12.79);
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
    ppi_post_cfar_threshold(17.5, 2e-20, 12.79);
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
          timing2.contains("last 200 of 1200 PRIs") && timing2.contains("not part of Total"),
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
        auto noise = [&seed]() { seed = seed * 1664525u + 1013904223u; return 800.0f * (0.8f + 0.4f * (seed >> 8) / 16777216.0f); };
        for (double az = 30.0; az < 70.0; az += 0.036) {
            for (int i = 0; i < nBins; i++) row[i] = i < 55 ? 0.0f : noise();   /* no data inside 4.1 km */
            if (az > 45.0 && az < 55.0)
                for (int i = 400; i < 404; i++) row[i] = 40000.0f;
            ppi_post_video(row.data(), nBins, dR, az, 20.0 + az / 36.0);
        }
        ppi_post_antenna(70.0, 22.0, 2);
        win.refresh();
        check(pw.videoBinsFilled() == 160 && std::fabs(pw.videoNoise_dB() - 10.0 * std::log10(672.0)) < 0.3,
              "video of 30 ... 70 deg (160 azimuth bins), noise floor (10th percentile) of the rows");
        const double r = 30.1e3;
        check(pw.videoIndexAt(r * std::sin(50 * kDeg), r * std::cos(50 * kDeg)) > 80 &&
              pw.videoIndexAt(20e3 * std::sin(40 * kDeg), 20e3 * std::cos(40 * kDeg)) < 20 &&
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
            ppi_post_video(row.data(), nBins, dR, az + 180.0, 30.0 + az / 36.0);           /* video: dropped */
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
        check(win.tabs()->count() == 2 && win.tabs()->currentIndex() == 0 && ppi_rdmap_wanted(&bm) == 0 &&
              pw.parentWidget() == win.tabs()->widget(0), "two tabs; the map is computed only in the second");
        win.resize(1400, 800);
        win.tabs()->setCurrentIndex(1);
        QCoreApplication::processEvents();
        check(ppi_rdmap_wanted(&bm) == 1 && bm == 1 && pw.parentWidget() == win.tabs()->widget(1) &&
              !win.rd().hasFrame(), "second tab: map of beam 2 wanted, the PPI moved next to it");
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
        m.beam = 2;  m.time_s = 41.0;
        ppi_post_rdmap(&m);                                       /* a map of another beam is not shown */
        win.refresh();
        check(rdw.frame().time_s == 40.0, "a map of another beam is not taken");
        if (!folder.empty()) {
            const QString shot = QString::fromStdString(folder + "/ppi_rd_tab.png");
            check(win.saveSnapshot(shot) && QImage(shot).width() == win.tabs()->currentWidget()->width(),
                  "screenshot of the second tab: map and small PPI");
        }
        win.tabs()->setCurrentIndex(0);
        check(ppi_rdmap_wanted(nullptr) == 0 && pw.parentWidget() == win.tabs()->widget(0), "back to the PPI tab");
    }

    if (!folder.empty()) {
        const QString shot = QString::fromStdString(folder + "/ppi_window_screenshot.png");
        check(win.ppi().saveSnapshot(shot) && QImage(shot).width() == win.ppi().width(),
              "screenshot of the PPI (without the panel)");
        QCoreApplication::processEvents();               /* layout after the text changes */
        win.grab().save(QString::fromStdString(folder + "/ppi_window.png"));
    }
    std::printf(failures ? "%d FAILURE(S)\n" : "PASS\n", failures);
    return failures ? 1 : 0;
}

use eframe::egui;
use hmc::gui::HmcApp;

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1280.0, 800.0])
            .with_min_inner_size([800.0, 500.0])
            .with_title("HMC Sandbox"),
        ..Default::default()
    };

    eframe::run_native(
        "hmc_gui",
        options,
        Box::new(|cc| Ok(Box::new(HmcApp::new(cc)))),
    )
}

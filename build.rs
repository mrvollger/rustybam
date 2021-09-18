// in build.rs
use clap::{crate_version, load_yaml, App, AppSettings};
use clap_generate::{
    generate, generate_to,
    generators::{Bash, Zsh},
};


fn main() {
    let yaml = load_yaml!("src/cli.yaml");
    let mut app = App::from(yaml)
        .version(crate_version!())
        .setting(AppSettings::SubcommandRequiredElseHelp);

    app.set_bin_name("rustybam");
    let outdir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("completions/");

    generate_to::<Bash, _, _>(&mut app, "rustybam", &outdir)
        .expect("Failed to generate bash completions");
    generate_to::<Zsh, _, _>(&mut app, "rustybam", &outdir)
        .expect("Failed to generate zsh completions");
}

// in build.rs
use clap_generate::{generate_to, generators::*};
include!("src/cli.rs");

fn main() {
    let mut app = Cli::into_app();
    let outdir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target/");
    /*
    // creates an error with cargo publish and then installing
    for bin_name in ["rb", "rustybam"] {
        // Bash, Fish, Zsh, PowerShell, Elvish
        generate_to(Bash, &mut app, bin_name, &outdir).unwrap();
        generate_to(Fish, &mut app, bin_name, &outdir).unwrap();
        generate_to(Zsh, &mut app, bin_name, &outdir).unwrap();
        generate_to(PowerShell, &mut app, bin_name, &outdir).unwrap();
        generate_to(Elvish, &mut app, bin_name, &outdir).unwrap();
    }
    */
}

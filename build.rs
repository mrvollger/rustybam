// in build.rs
use clap::{crate_version, load_yaml, App, AppSettings};
use clap_generate::{
    generate,
    generators::{Bash, Zsh},
};
use std::io;

fn main() {
    let yaml = load_yaml!("src/cli.yaml");
    let mut app = App::from(yaml)
        .version(crate_version!())
        .setting(AppSettings::SubcommandRequiredElseHelp);

    generate::<Bash, _>(&mut app, "rustybam", &mut io::stdout());
    generate::<Zsh, _>(&mut app, "rustybam", &mut io::stdout());

    /*for man in clap_generate::gen gen_manuals(&app) {
        let name = "rustybam.1";
        let mut out = std::fs::File::create(name).unwrap();
        use std::io::Write;
        out.write_all(man.render().as_bytes()).unwrap();
    }*/
}

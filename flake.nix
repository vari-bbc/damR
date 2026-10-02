{
  description = "Flake to get metaboapp development environment";
  inputs.nixpkgs.url = "github:rstats-on-nix/nixpkgs/r-bioc-devel";
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = nixpkgs.legacyPackages.${system};

      Imports = with pkgs.rPackages; [
        collapse
        ComplexHeatmap
        data_table
        ggplot2
        limma
        plotly
        RColorBrewer
        rlang
        stringi
        tmvtnorm
        viridis
      ];

      Suggests = with pkgs.rPackages; [
        testthat
        patchwork
        ggiraph
      ];

      rDevDeps = with pkgs.rPackages; [
        covr
        devtools
        DT
        htmltools
        pkgdown
        roxygen2
        V8
      ];

      devDeps = with pkgs; [
        checkbashisms
        html-tidy
      ];

      # default package
      rBuildInputs = [ ];
      rPropagatedBuildInputs = Imports ++ Suggests;
      rNativeBuildInputs = with pkgs; [
        pkg-config
        openssl
      ];

      metaboapp = pkgs.rPackages.buildRPackage {
        pname = "metaboapp";
        version = "0.0.9000";
        src = self;
        nativeBuildInputs = rNativeBuildInputs;
        propagatedBuildInputs = rPropagatedBuildInputs;
      };

      # Create R development environment with iscream and other useful libraries
    in {
      packages.default = metaboapp;
      devShells.default = pkgs.mkShell {
          inputsFrom = pkgs.lib.singleton metaboapp;
          packages = pkgs.lib.singleton devDeps ++ rDevDeps;
          shellHook = ''
            mkdir -p "$HOME/.R"
            export R_LIBS_USER="$HOME/.R"
          '';
      };
    });
}

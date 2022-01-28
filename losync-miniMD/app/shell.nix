let
  pkgs = import (builtins.fetchTarball
    "https://pm.bsc.es/gitlab/rarias/bscpkgs/-/archive/nosv/bscpkgs-master.tar.gz");
in
  pkgs.mkShell {
    name = "minimd";
    buildInputs = with pkgs.bsc; [ nanos6 extrae impi icc mcxx tampi ];
    shellHook = ''
      export LANG=C
    '';
  }

with import <nixpkgs> {};
    stdenv.mkDerivation {
        name = "MPI_LeapFrop_Simulation";
        buildInputs = [openmpi gcc];
        src = ./.;

        buildPhase = ''
            mpicc main.c -Ofast -lm -ffast-math -march=native -o simulation  

        '';

        installPhase = ''
            
            mkdir -p $out/bin
            cp simulation $out/bin
        '';

        unpackPhase = "";
    }
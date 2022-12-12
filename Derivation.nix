with import <nixpkgs> {};
    stdenv.mkDerivation {
        name = "MPI_LeapFrop_Simulation";
        buildInputs = [openmpi gcc];
        src = ./.;

        buildPhase = ''
            mpicc main.c -O3 -lm -o simulation  

        '';

        installPhase = ''
            
            mkdir -p $out/bin
            cp simulation $out/bin
        '';

        unpackPhase = "";
    }
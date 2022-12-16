{pkgs ? import <nixpkgs> {}}:
with pkgs;
let
  my-python-packages = python-packages: with python-packages; [
    pandas
    matplotlib
    numpy
    # other python packages you want
  ]; 
  python-with-my-packages = python3.withPackages my-python-packages;
in 
pkgs.mkShell {
    nativeBuildInputs = [pkgs.ccls pkgs.openmpi pkgs.gcc python-with-my-packages];
}
{ pkgs ? import (builtins.fetchTarball {
    url = "https://github.com/NixOS/nixpkgs/archive/nixos-unstable.tar.gz";
  }) {} }:

pkgs.mkShell {
  packages = with pkgs; [
    gcc
    gfortran
    gnumake
    openmpi
    cmake
    python3
  ];
}

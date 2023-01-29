#!/usr/bin/python
from subprocess import call
from os.path import expanduser
from os import chdir

def setupMachine():
  home_folder = expanduser("~")

  # Install useful programs
  command = ["sudo", "apt-get", "install"]
  command.extend(["git", \
                  "tig", \
                  "g++", \
                  "libgflags-doc", \
                  "libgflags-dev", \
                  "libgflags2", \
                  "make", \
                  "cmake", \
                  "python-pip", \
                  "htop", \
                  "openjdk-7-jdk", \
                  "gradle", \
                  "autoconf2.13", \
                  "zlib1g-dev", \
                  "libsasl2-dev", \
                  "libssl-dev", \
                  "gdb", \
                  "libhdf5-serial-dev", \
                  "libcr-dev", \
                  "mpich2", \
                  "mpich2-doc"
                  ])
  call(command)

def main():
  setupMachine()

if __name__ == "__main__":
  main()
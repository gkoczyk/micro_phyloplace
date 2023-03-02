#! /bin/bash
docker run -v `pwd`:/current  -e CUSERID=$(id -u) -e CGROUPID=$(id -g) -t --rm --detach=false micro_phyloplace $@

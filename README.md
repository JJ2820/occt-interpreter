# occt-interpreter

occt version - https://git.dev.opencascade.org/gitweb/?p=occt.git;a=snapshot;h=80ffc5f84dae96de6ed093d3e5d2466a9e368b27;sf=tgz

# docker commands

docker build -t xibyte/occt:wasm-builder_1.0 .

docker run -it -v $(pwd)/build-wasm:/build -v $(pwd)/occt:/occt xibyte/occt:wasm-builder_1.0 -i

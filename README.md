
Build and Run GUI:
-------
You need [CMake], [Eigen]v3, and [GTK+]v2.

    git clone https://github.com/aashi7/chomp_edt.git
    mkdir chomp_edt/build
    cd chomp_edt/build
    cmake ..
    make
    ./demo

The obstacle cost function is based on Euclidean distance transform. Environment thus needs to be discretized.
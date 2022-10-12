
**Berglund**:

[Installation]:

Dependencies: 
Boost, version 1.58.0 or newer
Armadillo, version 6.400 or newer

Build and install the code as follows. Note that `<INSTALL_PREFIX>` is a path below which the program will be installed. This could be e.g. `$HOME/local` to install into a user-local prefix.
```{c}
cd build
./gen_build.sh -DCMAKE_INSTALL_PREFIX=<INSTALL_PREFIX>
make
make install
```

The above will build both a release and a debug version of the code. Please use `make release` or `make debug` in place of `make` above if you want to build only the release or debug version. The binary for the release version will be called `std` and the binary for the debug version will be called `std-dbg`.

Note that `<INSTALL_PREFIX>/bin` and `<INSTALL_PREFIX>/lib` have to be included in your `PATH` and `LD_LIBRARY_PATH` environment variables, respectively. You can use 
```{c}
export PATH=<INSTALL_PREFIX>/bin:$PATH
export LD_LIBRARY_PATH=<INSTALL_PREFIX>/lib:$LD_LIBRARY_PATH
```


[Run]:
We can run the binary `std` on our seqFISH data.
```{c}
./std --file ../datasets/seqFISH+/counts.txt --output ../seqFISH_10000_Result/
```
Note that the inputs and outputs filenames here may require minor changes 
if you want to run your own examples. Please see the documents of the 
binary to get an idea of what a input file would look like. 
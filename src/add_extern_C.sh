#!/bin/bash

# This script modifies headers to allow C++ linking.

# cd rolloffdir
for header_file in *.h nicksrc/*.h
do
    echo $header_file
    mv $header_file $header_file.bak
    echo '#ifdef __cplusplus
extern "C" {
#endif' | cat - $header_file.bak > $header_file
    echo '#ifdef __cplusplus
}
#endif' >> $header_file
done

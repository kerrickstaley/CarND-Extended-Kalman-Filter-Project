#!/bin/bash
clang++ tools.cpp test_tools.cpp /usr/src/gtest/src/gtest-all.cc -I/usr/src/gtest -lpthread -o test_tools
./test_tools

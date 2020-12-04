#!/bin/bash

rivet-mkhtml --pwd Rivet_pp200.yoda
cp -r rivet-plots/ rivet-plots-pp200/

rivet-mkhtml --pwd Rivet_pp62.yoda
cp -r rivet-plots/ rivet-plots-pp62/


#!/usr/bin/env bash

grids=500
degree_of_polygons=(6)
nDs=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1)
n0s=(100)
Hcs=(-1)
vector_of_side_lengths=("[4.5 3]")

inputFolder="./inputFiles/"
inputFile="input_parameters.m"
rm -rf $inputFolder
mkdir -p $inputFolder

num=0
# non-exchange
for degree_of_polygon in ${degree_of_polygons[@]}; do
    for nD in ${nDs[@]}; do
        for n0 in ${n0s[@]}; do
            for Hc in ${Hcs[@]}; do
                for side_length in ${!vector_of_side_lengths[@]}; do
                    [[ ${vector_of_side_lengths[${side_length}]} =~ "[(*)\ (*)]" ]] && echo ${BASH_REMATCH}
#                    vector_1=${BASH_REMATCH[0]}
#                    vector_2=${BASH_REMATCH[1]}
                    sed -e "/^degree_of_polygon/c\ 
degree_of_polygon = ${degree_of_polygon};" \
                        -e "/^n_D/c\ 
n_D = ${nD};" \
                        -e "/^n0/c\ 
n0 = ${n0};" \
                        -e "/^critical_height/c\ 
critical_height = ${Hc};" \
                        -e "/^vector_of_side_lengths/c\ 
vector_of_side_lengths = ${vector_of_side_lengths[${side_length}]};" \
                        -e "/^number_of_triangles/c\ 
number_of_triangles = ${grids};" \
                        -e "/^include_ex/c\ 
include_ex = false;" \
                        -e "/^draw_ex/c\ 
draw_ex = false;" \
                        ${inputFile} > "${inputFolder}input_parameters_${num}.m"
                    num=$(( ${num}+1 ))
                    echo ${num}
                done
            done
        done
    done
done


# exchange
for degree_of_polygon in ${degree_of_polygons[@]}; do
    for nD in ${nDs[@]}; do
        for n0 in ${n0s[@]}; do
            for Hc in ${Hcs[@]}; do
                for side_length in ${!vector_of_side_lengths[@]}; do
                    [[ ${vector_of_side_lengths[${side_length}]} =~ "[(*)\ (*)]" ]] && echo ${BASH_REMATCH}
#                    vector_1=${BASH_REMATCH[0]}
#                    vector_2=${BASH_REMATCH[1]}
                    sed -e "/^degree_of_polygon/c\ 
degree_of_polygon = ${degree_of_polygon};" \
                        -e "/^n_D/c\ 
n_D = ${nD};" \
                        -e "/^n0/c\ 
n0 = ${n0};" \
                        -e "/^critical_height/c\ 
critical_height = ${Hc};" \
                        -e "/^vector_of_side_lengths/c\ 
vector_of_side_lengths = ${vector_of_side_lengths[${side_length}]};" \
                        -e "/^number_of_triangles/c\ 
number_of_triangles = ${grids};" \
                        -e "/^include_ex/c\ 
include_ex = true;" \
                        -e "/^draw_ex/c\ 
draw_ex = true;" \
                        ${inputFile} > "${inputFolder}input_parameters_${num}.m"
                    num=$(( ${num}+1 ))
                    echo ${num}
                done
            done
        done
    done
done

echo "Total" ${num} "files."

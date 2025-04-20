//
// Created by 阮立心 on 25-4-19.
//

#ifndef GRID_DATA_H
#define GRID_DATA_H
#include <vector>

#include "type.h"

using namespace CGL;

template <typename T>
class GridData {
public:
    const unsigned int nRows;
    const unsigned int nCols;

    std::vector<T> data;

    GridData(const unsigned int _nRows, const unsigned int _nCols)
        : nRows(_nRows), nCols(_nCols) {
        data.resize(nRows * nCols);
    }

    T get(unsigned int i, unsigned int j) const { return data[j + nCols * i]; }

    void set(unsigned int i, unsigned int j, const T& value) {
        data[j + nCols * i] = value;
    }
};
#endif //GRID_DATA_H

//
// Created by martin on 10/12/2019.
//

#ifndef SWE_FLOAT2D_BUFFER_HPP
#define SWE_FLOAT2D_BUFFER_HPP

#include <memory>

#include "types/Float2D.hpp"
#include "types/Float2DNative.hpp"

class Float2DBuffer : public Float2D {
public:
    Float2DBuffer() :
            Float2D(0, 0) {};

    Float2DBuffer(int cols, int rows, bool localTimestepping, Float2DNative &realData) :
            Float2D(cols, rows) {

        if (localTimestepping) {

            std::shared_ptr<float> tmp(new float[rows * cols], std::default_delete<float[]>());
            data = tmp;
            rawData = data.get();

        } else {
            // If there is no local timestepping buffer points to h |hu | hv
            data = realData.getPointer();
            rawData = realData.getPointer().get();
        }

    }

    ~Float2DBuffer() {}

private:


    std::shared_ptr<float> data;
};

#endif // SWE_FLOAT2D_BUFFER_HPP

#ifndef SWE_FLOAT2D_NATIVE_HPP
#define SWE_FLOAT2D_NATIVE_HPP

#include <memory>

#include "types/Float2D.hpp"

class Float2DNative : public Float2D
{
public:
    Float2DNative() : Float2D(0, 0){};

    Float2DNative(int cols, int rows) : Float2D(cols, rows)
    {
        std::shared_ptr<float> tmp(new float[rows * cols], std::default_delete<float[]>());
        data    = tmp;
        rawData = data.get();
    }

    Float2DNative(const Float2DNative& src) : Float2D(src.cols, src.rows)
    {
        std::shared_ptr<float> tmp(new float[rows * cols], std::default_delete<float[]>());
        data    = tmp;
        rawData = data.get();
        memcpy(rawData, src.rawData, cols * rows);
    }

    ~Float2DNative() {}

    std::shared_ptr<float> getPointer()
    {
        return data;
    }

private:
    std::shared_ptr<float> data;
};

#endif // SWE_FLOAT2D_NATIVE_HPP

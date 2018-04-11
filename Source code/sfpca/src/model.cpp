#include "model.h"

MODEL_TYPE string_to_MT(std::string model_type_string)
{
    if (model_type_string.compare("PCA") == 0)
        return PCA;
    else if (model_type_string.compare("LDA") == 0)
        return LDA;
    else if (model_type_string.compare("PLS") == 0)
        return PLS;
    else if (model_type_string.compare("CCA") == 0)
        return CCA;
    else
        throw std::invalid_argument(model_type_string + " is not currently supported");
}

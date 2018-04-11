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

Model build_model(arma::mat X, arma::mat Y, std::string model_type)
{
    Model model = {.model_type = string_to_MT(model_type)};

    switch (model.model_type)
    {
    case PCA:
        if (DEBUG)
            cout << "Chosing case PCA: " << model.model_type << endl;

        model.X = X;
        break;
    case LDA:
        throw std::invalid_argument(model_type + " is not currently supported");
        break;
    case PLS:
        model.X = X.t() * Y;
        break;
    case CCA:
        model.X = X.t() * Y;
        break;
    default:
        throw std::invalid_argument(model_type + " is not currently supported");
    }
    return model;
};
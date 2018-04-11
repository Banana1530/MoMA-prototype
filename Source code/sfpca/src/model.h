#ifndef _MODEL_H_
#define _MODEL_H_
#include "common.h"
typedef enum { PCA,
               LDA,
               CCA,
               PLS } MODEL_TYPE;
typedef struct
{
    MODEL_TYPE model_type;
    arma::mat X;
    bool is_sparse;
} Model;

MODEL_TYPE string_to_MT(std::string model_type_string);
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
#endif  
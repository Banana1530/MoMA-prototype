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
Model build_model(arma::mat X, arma::mat Y, std::string model_type);
#endif  
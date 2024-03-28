#include <iostream>
#include "legendre.h"
#include <vector>

int main(){
    std::vector<int> v{ 1, 5, 6, 7, 1, 3, 4, 5, 3, 6};
    for ( auto &c : v){
        std::cout<<c<<"\t";
    }

    

    std::cout<<"\n--------\n";
    for ( size_t i = 0; i <v.size(); i++ ){
        //std::cout<<"--> "<<v[i]<<"\n";
            for (size_t j=i; j < v.size(); j++){

                if ( v[i] == v[j] && (i!=j) ){
                    std::cout<<" i = "<<i<<", j = "<<j<<", value = "<<v[i]<<"\n";
                }
            }
        }


}
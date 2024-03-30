#include <iostream>
#include <vector>
#include </home/singhp10/MonteCarlo/abeille/build/_deps/yaml-cpp-src/include/yaml-cpp/yaml.h>

int main(){
    YAML::Node input_ = YAML::LoadFile("file.yaml");

    if ( input_["person"]) {
        for(size_t it = 0; it < input_["person"].size(); it++){
            std::string name = input_["person"][it]["name"].as<std::string>();
            int age = input_["person"][it]["age"].as<int>();
            std::cout<<"Name "<< name<<", age = "<<age<<"\n";
        } 
    }else{
        std::cout<<"Not found"<<"\n";
    }

    std::vector<int> ve = input_["list"].as<std::vector<int>>();
    for ( auto &c: ve){
        std::cout<<c<<"\t";
    }
    std::cout<<"\n";
}

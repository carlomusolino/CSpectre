#include <vector>
#include <iostream>

using std::cout ;

int main() {

    std::vector<double> x;

    cout << x.size() << "\n";

    for(int i=0; i<10; i++){
        x.push_back(2.0*i);
    }

    for (const auto& val : x) cout << val << "\n";
    
    

}
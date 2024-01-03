#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <bitset>
#include <iomanip>
#include <time.h>
#include <chrono>
#include <cstdint>
using namespace std;

class MBitInteger {
private:
    vector<uint64_t> data;
    bool Isnegative = false;

public:

    MBitInteger(){};
    
    
    MBitInteger(const string& str) { 
        string bit;
        uint64_t intbit;
        int i = 16;
        // cout << str.length() << endl;
        while(i < str.length()){
            bit = str.substr(str.length()-i,16);
            stringstream(bit) >> hex >> intbit;
            data.insert(data.begin(),intbit);
            // cout << "bot bit = " << bit << endl;
            i=i+16;
        }
        
        bit = str.substr(0,str.length()-(i-16));
        // cout << "bot bit2 = " << bit << endl;
        stringstream(bit) >> hex >> intbit;
        data.insert(data.begin(),intbit);

        // cout << "--------our result----------" << endl;
        //     for(int i = 0; i < data.size();i++){
        //         cout << data[i] << endl;
        //     }
        //     cout << "--------end-----------" << endl;

    }
    
   
    
    void print(){
        // cout << "OUR NUBER IS: " << endl;
        // for(int i = 0; i < data.size();i++){
        //     cout << data[i] << endl;
        // }
        // cout << "in 16: " << endl;
        cout << hex << data[0];
        for(int i = 1; i < data.size(); i++){
            cout << hex << setw(16) << setfill('0') << data[i];
        }
        cout << dec;
        cout << endl;
    }
    
    
    
    int maxn(uint64_t number){
        uint64_t result = 1;
        int i = 0;
        
        while(result <= number && i != 20){
            
            result *= 10;
            i++;
        }
        cout << "number " << number << " and i " << i << endl;
        return i;
    }
    
    int symb(uint64_t number,int order){
        int result;
        int i = 0;
        while(i != order-1 && number != 0){
            number = number/10;
            i++;
        }
        
        if(number == 0){
            
            return 0;
        }else{
            
            return number%10;
        }
        
    }
    
    
    
    void interpret(int sz){
        int dif = sz-data.size();
        if(dif > 0){
            data.resize(sz,0);
            rotate(data.begin(), data.begin()+(sz-dif), data.end());
        }else if(dif < 0){
            data.resize(sz);
        }
    }
    void interpret(){
        int dif = 0;
        for(int i = 0; i < data.size();i++){
            if(data[i] != 0 ||(i == data.size()-1 && data[i] == 0)){
                rotate(data.begin(), data.begin()+dif, data.end());
                data.resize(data.size()-dif);
                break;
            }else{
                dif++;
            }
        }
        
    }
    
    bool operator>=(MBitInteger& second){
        if(data.size() > second.data.size()){
            return true;
        }else if(data.size() < second.data.size()){
            return false;
        }
        for (int i = 0; i < data.size();i++) {
            if (data[i] < second.data[i]) {
                return false;
            }else if (data[i] > second.data[i]) {
                return true;
            }
        }
        return true;
    }
    bool operator>(MBitInteger& second){
        if(data.size() > second.data.size()){
            return true;
        }else if(data.size() < second.data.size()){
            return false;
        }
        for (int i = 0; i < data.size();i++) {
            if (data[i] < second.data[i]) {
                return false;
            }else if (data[i] > second.data[i]) {
                return true;
            }
        }
            return false;
    }
    
    bool operator<=(MBitInteger& second){
        if(data.size() < second.data.size()){
            return true;
        }else if(data.size() > second.data.size()){
            return false;
        }
        for (int i = 0; i < data.size();i++) {
            if (data[i] > second.data[i]) {
                return false;
            }else if (data[i] < second.data[i]) {
                return true;
            }
        }
        return true;
    }
    bool operator<(MBitInteger& second){
        if(data.size() < second.data.size()){
            return true;
        }else if(data.size() > second.data.size()){
            return false;
        }
        for (int i = 0; i < data.size();i++) {
            if (data[i] > second.data[i]) {
                return false;
            }else if (data[i] < second.data[i]) {
                return true;
            }
        }
        return false;
    }
    
    
    
    MBitInteger operator+(MBitInteger& second) {
        MBitInteger result;
        int maxm = max(data.size(),second.data.size());
        result.data.resize(maxm+1,0);
        uint64_t carry = 0;
        interpret(maxm);
        second.interpret(maxm);
        for(int i = maxm-1; i >= 0;i--){
            result.data[i] = (data[i] + second.data[i] + carry);
            if(result.data[i] < data[i] || result.data[i] < second.data[i]){
                carry = 1;
            }else{
                carry = 0;
            }
        }
        if(carry == 1){
            rotate(result.data.begin(), result.data.begin()+maxm, result.data.end());
            result.data[0] = carry;
        }else{
            result.data.resize(maxm);
        }
        
        return result;
    }
    MBitInteger operator-(MBitInteger& second) {
        MBitInteger result;
        int maxm = max(data.size(),second.data.size());
        result.data.resize(maxm,0);
        interpret(maxm);
        second.interpret(maxm);
        MBitInteger swap;
        MBitInteger swap2;
        
        for(int i = 0;i < maxm;i++){
            if(data[i] < second.data[i]){
                swap.data = second.data;
                swap2.data = data;
                result.Isnegative = true;
                break;
            }else if(data[i] > second.data[i]){
                swap2.data = second.data;
                swap.data = data;
                break;
            }else if(i == maxm-1 && data[i] == second.data[i]){
                result.interpret();
                return result;
            }
        }

        uint64_t borrow = 0;
        int temp = 0;
        for(int i = maxm-1; i >= 0;i--){
            result.data[i] = swap.data[i] - swap2.data[i] - borrow;
            if(swap.data[i] < swap2.data[i]+borrow){
                borrow = 1;
            }else{
                borrow = 0;
            }
        }
        
        return result;
    }
    
    MBitInteger operator*(MBitInteger& second){
        // cout << "one*two: ";
        // cout << "---------------" << endl;
        // print();
        // second.print();
        MBitInteger result;
        result.data.resize(data.size()+second.data.size());
        uint64_t carry = 0;
        uint64_t carry2 = 0;
        uint64_t ost[2];
        uint64_t num[2];
        uint64_t mult[4];
        uint64_t ostmult[2] = {0,0};
        uint64_t u32 = 4294967296;
        uint64_t help1;
        uint64_t help2;
        uint64_t before;
        
        for(int i = data.size()-1;i >= 0;i--){
            for(int j = second.data.size()-1;j >= 0;j--){
                if(data[i] < u32 && second.data[j] < u32){
                    // cout << "a123b123c123" << endl;
                    // cout << i << " " << j << " " << data[i]*second.data[j] << " "  << result.data[i+j+1] << " " << carry << " " << result.data[0] << " " << carry2 << endl;
                    result.data[i+j+1] += data[i]*second.data[j];
                    // result.data[i+j+1] += carry;
                    if(result.data[0] > carry || result.data[0] > carry2){
                        carry = 1;
                    }else{
                        carry = 0;
                    }
                    
                    result.data[i+j+1] += carry;
                    break;
                }
                if(data[i] >= u32){
                    num[0] = data[i]/u32;
                    ost[0] = data[i]%u32;
                    help1 = 1;
                }else{
                    num[0] = data[i];
                    ost[0] = 0;
                    help1 = u32;
                }
                if(second.data[j] >= u32){
                    num[1] = second.data[j]/u32;
                    ost[1] = second.data[j]%u32;
                    help2 = 1;
                }else{
                    num[1] = second.data[j];
                    ost[1] = 0;
                    help2 = u32;
                }
                
                mult[0] = num[0]*num[1];
                mult[1] = num[0]*ost[1]/u32;
                ostmult[0] = (num[0]*ost[1])%u32;
                mult[2] = num[1]*ost[0]/u32;
                ostmult[1] = (num[1]*ost[0])%u32;
                mult[3] = ((ost[0]*ost[1]/u32)+ostmult[0]+ostmult[1])/u32;
                before = result.data[i+j+1];
                result.data[i+j+1] += data[i]*second.data[j];
                carry = 0;
                for(int x = 0;x < 4;x++){
                    carry += mult[x];
                }
                carry = carry/help1/help2;
                
                
                if(result.data[i+j+1] < before){
                    result.data[i+j] += 1;
                }
                before = result.data[i+j];
                result.data[i+j] += carry;
                if(result.data[i+j] < before){
                    result.data[i+j-1] += 1;
                }
                
            }
            
        }
        // cout << "-------------------------" << endl;
        // cout << result.data[5] << endl;
        result.interpret();
        return result;
    }
    
    
    
    MBitInteger karatsuba(MBitInteger& x, MBitInteger& y,int num) {
    size_t n = max(x.data.size(), y.data.size());
    
    x.interpret(n);
    y.interpret(n);
        // Базовий випадок: якщо числа коротші за деяку константу, використовувати звичайне множення
        if (n <= 1) {
            
            return x * y;
            
        }

        size_t half_size = n / 2;
    
        // Розбиваємо числа на дві половини
        MBitInteger xl = x.lowHalf(half_size);
        MBitInteger xr = x.highHalf(half_size);
        MBitInteger yl = y.lowHalf(half_size);
        MBitInteger yr = y.highHalf(half_size);
        // cout << "-----------" << endl;
        // cout << "num: " << num << endl;
        // cout << "-----------" << endl;
        // xl.print();
        // xr.print();
        // cout << endl;
        // yl.print();
        // yr.print();
        // Рекурсивно обчислюємо три проміжні добутки
        MBitInteger P1 = karatsuba(xl, yl,num+1);
        MBitInteger P2 = karatsuba(xr, yr,num+1);
        MBitInteger sumx = xl + xr;
        MBitInteger sumy = yl + yr;
        MBitInteger P3 = karatsuba(sumx, sumy,num+1);
        // cout << "--------------------------------------------" << endl;
        // xl.print();
        // xr.print();
        // cout << endl;
        // yl.print();
        // yr.print();
        // cout << "sumx: ";
        // sumx.print();
        // cout << "sumy: ";
        // sumy.print();
        // cout << "P1: ";
        // P1.print();
        // cout << "P2: ";
        // P2.print();
        // cout << "P3: ";
        // P3.print();
        
        // Обчислюємо кінцевий результат
        // MBitInteger P4;
        // P4.data.resize(n);
        // P4.data[0] = 1;
        
        // cout << "P4: ";
        // P4.print();
        // MBitInteger P5;
        // P5.data.resize((n)/2);
        // P5.data[0] = 1;
        // cout << "P5: ";
        // P5.print();
        // cout << n << " size" << endl;
        // cout << "P4: ";
        // P4.print();
        // cout << "P1: ";
        // P1.print();
        MBitInteger result = P1;
        if(n%2==1 ){
            n++;
        }
        result.data.resize(result.data.size()+n);
        // cout << "result: ";
        // result.print();
        MBitInteger P = P3 - P1;
        P = P - P2;
        
        // cout << "P ";
        
        // P.print();
        // P = P5 * P;
        if(n%2==1 ){
            n++;
        }
        P.data.resize(P.data.size()+n/2);
        // P.print();
        // result = result + ((P3 - P1 - P2) * (1ULL << half_size));
        result = result + P;
        // cout << "result+P: ";
        // result.print();
        // result.data.resize(result.data.size()+n-1);
        // result = P4 * result ;
        // cout << "result*P4: ";
        // result.print();
        result = result + P2;
        // cout << "result+P2: ";
        // result.print();
        
        // cout << "ответ ";
        // result.print();
        // P1.print();
        // P.print();
        // P2.print();
        // cout << "____________________" << endl;
        return result;
}
    
    
    MBitInteger highHalf(size_t half_size) const {
        MBitInteger result(*this);
        result.data.erase(result.data.begin(), result.data.begin() + min(result.data.size(), half_size));
        return result;
    }

    MBitInteger lowHalf(size_t half_size) const {
        MBitInteger result(*this);
        result.data.erase(result.data.begin() + min(result.data.size(), half_size), result.data.end());
        return result;
    }
    
    
    
    
    
    MBitInteger operator/(MBitInteger& second) {
        vector<uint64_t> B;
        vector<uint64_t> A;
        MBitInteger Q;
        Q.data.resize(data.size());
        A = bit2();
        B = second.bit2();
        vector<uint64_t> R;
        vector<uint64_t> C;
        int k = B.size();
        R = A;
        
        MBitInteger result;
        MBitInteger plus;
        plus.data.push_back(0);
        result.data.resize(A.size());
        
        int i = 0;

        while (!bit2compare(R,B)) {
            // cout << "our R: " << endl;
            // R.print();
            // cout << "------------" << endl << "our A: " << endl;
            // print();
            // cout << "------------" << endl << "our B: " << endl;
            // second.print();
            int t = R.size();
            C = B;
            C.resize(R.size());
            
            // cout << endl << "our C:" << endl;
            // for(int i = 0; i < C.size();i++){
            //     cout << C[i] << " ";
            // }
            
            // cout << endl << "our C2:" << endl;
            // for(int i = 0; i < C.size();i++){
            //     cout << C[i] << " ";
            // }
            // cout << endl;
            
            
            if (bit2compare(R,C)) {
                t = t - 1;
                rotate(C.begin(),C.begin()+C.size()-1,C.end());
                // cout << endl << "our C3:" << endl;
                // for(int i = 0; i < C.size();i++){
                //     cout << C[i] << " ";
                // }
                // C.resize(C.size()-1);
                // cout << endl << "our C3:" << endl;
                // for(int i = 0; i < C.size();i++){
                //     cout << C[i];
                // }
                // cout << endl;
                // cout << "------------" << endl << "R < C so our new C: " << endl;
                // C.print();
                
            }
            // cout << "our R before minus" << endl;
            // for(int i = 0; i < R.size();i++){
            //         cout << R[i];
            //     }
            // cout << endl;
            // cout << endl << "our C:" << endl;
            // for(int i = 0; i < C.size();i++){
            //     cout << C[i];
            // }
            // cout << endl;
            R = bitminus(R,C);
            // cout << "our R after minus" << endl;
            // for(int i = 0; i < R.size();i++){
            //         cout << R[i];
            //     }
            // cout << endl;
            
            // result.data = Q.data;
            i = (t-k)/64;
            // cout << "our result:" << endl;
            
            // cout << "t-k: " << t-k << endl;
            Q.data[Q.data.size()-i-1] += (1ULL << (t-k));
            // Q.print();
        }
        
        Q.interpret();
        
        return Q;
    }
    

    
    vector<uint64_t> bitminus(const vector<uint64_t>& R, const vector<uint64_t>& C) {
        if (R.size() != C.size()) {
            throw std::invalid_argument("Длины векторов должны быть одинаковыми");
        }
        
        std::vector<uint64_t> result(R.size(), 0);
        if(R >= C && C >= R){
            return result;
        }
        int borrow = 0;
        for (int i = R.size()-1; i >= 0; i--) {
            
            int diff = R[i] - C[i] - borrow;
            if (diff < 0) {
                diff += 2;
                borrow = 1;
            } else {
                borrow = 0;
            }
            result[i] = diff;
        }
        int i = 0;
        while(result[i] == 0){
            i++;
        }
        rotate(result.begin(),result.begin()+i,result.end());
        result.resize(result.size()-i);
        return result;
    }
    
    
    MBitInteger operator%(MBitInteger& second){
        MBitInteger ost,copy,dif;
        copy.data = data;
        MBitInteger result = copy/second;
        dif = second*result;
        ost = copy - dif;
        return ost;
    }
    
    MBitInteger operator^(MBitInteger& second){
        MBitInteger result;
        MBitInteger A;
        A.data = data;
        result.data = {1};
        vector<uint64_t> B = second.bit2();
        for(int i = 0; i < B.size(); i++){
            if(B[i] == 1){
                result = result*A;
            }
            A=A*A;
        }
        return result;
    }
    

    vector<uint64_t> bit2(){
        vector<uint64_t> result;
        for(int j = 0;j < data.size();j++){
            for (int i = sizeof(data[j]) * 8 - 1; i >= 0; --i) {
                int bit =  ((data[j] >> i) & 1);
                result.push_back(bit);
            }
        }
        while (!result.empty() && result.front() == 0) {
            result.erase(result.begin());
        }
        return result;
    }
    
    
    bool bit2compare(vector<uint64_t>& R, vector<uint64_t>& C){
        // cout << "our R:" << endl;
        // for(int i = 0; i < R.size();i++){
        //     cout << R[i];
        // }
        // cout << endl << "our C:" << endl;
        // for(int i = 0; i < C.size();i++){
        //     cout << C[i];
        // }
        // cout << endl;
        int i = 0;
        while(R[i] == 0 && i != R.size()){
            i++;
        }
        if(i == R.size()){
            return true;
        }
        if(R.size() == C.size()){
            
                for(int i = 0; i < R.size();i++){
                    
                    if(R[i] < C[i]){
                        return true;
                        
                    }else if(R[i] > C[i]){
                        return false;
                        
                    }
                }
                
                return false;
            }else if(R.size() < C.size()){
                return true;
            }else{
                return false;
            }
    }
    
    string vectorToString() const {
        stringstream ss;
        
        for (int i = 0; i < data.size(); i++) {
            ss << hex << setw(16) << setfill('0') << data[i];
        }
        return ss.str();
    }

    
    
    
    
    MBitInteger gcd(const MBitInteger& a, const MBitInteger& b) const {
        MBitInteger tempA = a;
        MBitInteger tempB = b;
        MBitInteger temp;
        MBitInteger ch;
        ch.data.push_back(0);

        while (tempB > ch) {
            temp.data = tempB.data;
            tempB = tempA % tempB;
            tempA.data = temp.data;
        }
        tempA.interpret();
        return tempA;
    }
    MBitInteger lcm(MBitInteger& a,MBitInteger& b){
        MBitInteger gcdAB = gcd(a, b);
        MBitInteger product = a * b;
        
        
        // НСК = (|a * b|) / НОД(a, b)
        return product / gcdAB;
    }
    

    
    MBitInteger mod(MBitInteger& a,MBitInteger& m) {
    // Вычисление параметров для алгоритма Барретта
    // MBitInteger mu("0");
    int k = (m.data.size()-1)*16;
    
    stringstream ss;
    ss << hex << setw(16) << m.data[0];
    string str = ss.str();
    str.erase(std::remove_if(str.begin(), str.end(), [](unsigned char c) { return std::isspace(c); }), str.end());
    // cout << "str " << str << " ss.str " << str.length() << " k " << k << endl;
    k += str.length();
    // k = 3;
    // k = std::ceil(static_cast<double>(k / 2.0));
    // int k = std::ceil(static_cast<double>(m.data.size()) / 2.0);
    // cout << "k " << k  << " m.data.size " << m.data.size() << endl;
    // cout << "a.data: ";
    // a.print();
    // cout << "m.data: ";
    // m.print();
    
    MBitInteger null;
    null.data.push_back(0);
    string smu = "1";
    for(int i = 0; i < 2*k;i++){
        smu += "0";
    }
    // cout << "smu " << smu << endl;
    MBitInteger mu(smu);
    // mu.data.resize(2*k);
    // mu.data[0] = 1;
    // cout << "mu.data: ";
    // mu.print();
    mu = mu / m;
    // cout << "mu2.data: ";
    // mu.print();
    
    // MBitInteger mu = (MBitInteger(1) << (2 * k)) / m;
    // MBitInteger mBar = (MBitInteger(1) << (2 * k)) % m;
    
    // Вычисление q
    MBitInteger q1 = a.KillLastDigits(k-1);
    // cout << "q1.data: ";
    // q1.print();
    q1 = q1 * mu;
    // cout << "q1.data: ";
    // q1.print();
    // MBitInteger q1 =https://www.onlinegdb.com/online_c++_compiler#tab-stdin (a >> (m.data.size() - k)) * mu;
    // MBitInteger q3 = a * q1;
    // cout << "q3.data: ";
    // q3.print();
    MBitInteger q2 = q1.KillLastDigits(k+1);
    // cout << "q2.data: ";
    // q2.print();
    // MBitInteger q2 = (a * q1) >> (m.size() + k);
    MBitInteger q = q2 * m;
    // cout << "q.data: ";
    // q.print();
    MBitInteger r = a - q;
    r.interpret();
    // m.print();
    
    
    // cout << "r.data: ";
    // r.print();
    
    int i = 0;
    while(r >= m && i < 2){
        r = r - m;
        i++;
    }
    m.interpret();
    // cout << i << endl;
    if(r >= m){
        // m.print();
        // cout << "-----------------------" << endl;
        r = mod(r,m);
    }
    

    return r;
}

    MBitInteger KillLastDigits(int count) const {
        MBitInteger result;
        string strRepresentation = vectorToString();
        for(int i = 0; i < count; i++){
            if (!strRepresentation.empty()) {
                strRepresentation.pop_back();
            }
        }
        // cout << "strRepresentation " << strRepresentation << endl;
        result = MBitInteger(strRepresentation);
        // result.interpret();
        // cout << "our result: ";
        // result.print();
        return result;
    }
    
    MBitInteger BarrettReduction(MBitInteger& A, MBitInteger& B, MBitInteger& N) {
        int k = (N.data.size()-1)*16;
        stringstream ss;
        ss << hex << setw(16) << N.data[0];
        string str = ss.str();
        str.erase(std::remove_if(str.begin(), str.end(), [](unsigned char c) { return std::isspace(c); }), str.end());
        // cout << "str " << str << " ss.str " << str.length() << " k " << k << endl;
        k += str.length();
        vector<uint64_t> B2 = B.bit2();
        MBitInteger C("1");
        MBitInteger CA;
        MBitInteger AA;
        string num = "1";
        for(int i = 0; i < 2*k; i++){
            num += "0";
        }
        MBitInteger u(num);
        u = u / N;
        // cout << "------------------" << endl;
        // for(int i = B2.size()-1;i >= 0; i--){
            
        //     cout << B2[i] << " ";
        // }
        // cout << endl << "------------------" << endl;
        for (int i = B2.size()-1;i >= 0; i--) {
            // cout << "===========================" << endl;
            if (B2[i] == 1) {
                CA = C*A;
                C = mod(CA, N);
            }
            AA = A * A;
            // cout << "===========================" << endl;
            A = mod(AA, N);
        }
    
        return C;
    }

    
    MBitInteger modplus(MBitInteger& one,MBitInteger& two,MBitInteger& m){
        MBitInteger result;
        result = one + two;
        result = mod(result,m);
        return result;
    }
    MBitInteger modminus(MBitInteger& one,MBitInteger& two,MBitInteger& m){
        MBitInteger result;
        result = one - two;
        result = mod(result,m);
        return result;
    }
    MBitInteger modmult(MBitInteger& one,MBitInteger& two,MBitInteger& m){
        MBitInteger result;
        result = one * two;
        result = mod(result,m);
        return result;
    }
    
    
};

int main() {
    double time_spent = 0.0;
    clock_t begin = clock();
    string str1 = "23C0D0050AE991EF232C32AA88639EF38290F68434288F9CBABABABABABABABABAABABE499B9220ABCA87E41472A1E15199B20E6CFDFE13CBB4E1F9358760F04A0FD9EA8A2B5A98E";
    string str2 = "1604077CEA1866B49C76BA21B38F3F0ABABABABABABABABABABAB6A4C5236EA0EB2B4C9BE1CF0CAD51CD5310E778A3A464E252B7215F3BE74DCFAB0C09DEDE3FA462A06C";
    // string str1 = "10000000000000000";
    // string str2 = "31280";
    
    // string str1 = "AAAA";
    // string str2 = "2";
    MBitInteger num1(str1);
    MBitInteger num2(str2);
    
    cout << "-------------MAIN CODE------------" << endl;
    cout << "Число 1: ";
    num1.print();
    cout << "Число 2: ";
    num2.print();
    
    cout << "Додавання ";
    MBitInteger result = num1 + num2;
    result.print();
    
    
    cout << "Різниця ";
    MBitInteger result2 = num1 - num2;
    result2.print();
    
    cout << "Множення ";
    auto startNaive = chrono::high_resolution_clock::now();
    MBitInteger result3 = num1 * num2;
    auto endNaive = std::chrono::high_resolution_clock::now();
    chrono::duration<double> durationNaive = endNaive - startNaive;
    result3.print();
    
    //cout << "4fe4a1a11946e17450e39277a652ad488fc5ab28e4c9c56bd4b338d327ce87af1d74e9f38853f18e848b0b3e9e6df6c5627445f34972472743489fe2d59f98d6ff879113c8217fb2ad6c13ec91a3051e8b74c8066a78db1092c03ec0846512fac85c84abe94960fcc77469853069501e0fb4ebb0a157d3760d5d03de0eb18cac4" << endl;
    cout << "Множення методом карацуби: ";
    auto startKaratsuba = std::chrono::high_resolution_clock::now();
    MBitInteger result30 = num1.karatsuba(num1,num2,1);
    auto endKaratsuba = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> durationKaratsuba = endKaratsuba - startKaratsuba;
    result30.interpret();
    result30.print();
    
    // string str3 = "12193263113";
    // MBitInteger num3(str3);
    // cout << "Число 3: ";
    // num3.print();
    
    // vector<uint64_t> res1 = num3.bit2();
    // cout << "Двійковий вигляд: ";
    // for(int i = 0;i< res1.size();i++){
    //     cout << res1[i];
    // }
    // cout << endl;
    // cout << "Ділення: ";
    // MBitInteger result4 = num1/num2;
    // result4.print();
    
    // cout << "Остача від Ділення: ";
    // MBitInteger ost1 = num1%num2;
    // ost1.print();

    // cout << "Степінь ";
    // MBitInteger result5 = num1*num1;
    // result5.print();
    
    // cout << "Перевірка на корректність "<< endl;
    // MBitInteger result6 = (num1 + num2) * num3;
    // result6.print();
    // MBitInteger result7 = (num1 * num3);
    // MBitInteger result8 = (num2 * num3);
    // MBitInteger result9 = result7 + result8;
    // result9.print();
    

    cout << "НСД чисел: ";
    MBitInteger result10 = num1.gcd(num1,num2);
    result10.print();
    cout << "НСК чисел: ";
    MBitInteger result11 = num1.lcm(num1,num2);
    result11.print();
    MBitInteger m("c202cd4c4f68a8d76ad1d");
    MBitInteger m2("c202cd4c4f68a8d76ad1d");
    cout << "m: ";
    m.print();
    // cout << "МОД^^^ чисел: ";
    // MBitInteger result15 = num1.mod(m2,m);
    // result15.print();
    cout << "МОД+ чисел: ";
    MBitInteger result12 = num1.modplus(num1,num2,m);
    result12.print();
    cout << "МОД^ чисел: ";
    MBitInteger result13 = num1.BarrettReduction(num1,num2,m);
    result13.print();
    cout << "МОД- чисел: ";
    MBitInteger result14 = num1.modminus(num1,num2,m);
    result14.print();
    cout << "МОД* чисел: ";
    
    MBitInteger result16 = num1.modmult(num1,num2,m);
    result16.print();
    
    
    
    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The elapsed time is %f seconds", time_spent);
    std::cout << "\nNaive multiplication Time: " << durationNaive.count() << " seconds" << std::endl;
    std::cout << "Karatsuba multiplication Time: " << durationKaratsuba.count() << " seconds" << std::endl;
    return 0;

    
    // cout << "Степінь ";
    // MBitInteger result5 = num1^num2;
    // result5.print();
    
    // cout << "Шіснадцядкове представлення: ";
    // MBitInteger result6 = num1.bit2();
    // result6.print();
    // // result6 = bit16(num2);
    // cout << result6 << endl;

    return 0;
    // cout << "Интерпретуємо!" << endl;
    // result.interpret(20);
    // result.print();
    // cout << "Вертаємо ";
    // result.interpret();
    // result.print();
    
}

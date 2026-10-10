#include <bits/stdc++.h>
using namespace std;

struct TermKey {
    int b, c, d; // exponents of x, sin, cos
    bool operator<(TermKey const& o) const {
        if (b != o.b) return b < o.b;
        if (c != o.c) return c < o.c;
        return d < o.d;
    }
    bool operator==(TermKey const& o) const {return b==o.b && c==o.c && d==o.d;}
};

using Poly = map<TermKey, long long>; // coefficient map

void poly_add(Poly &a, const Poly &b, long long mul=1) {
    for (auto &kv : b) {
        a[kv.first] += kv.second * mul;
        if (a[kv.first]==0) a.erase(kv.first);
    }
}

Poly poly_mul(const Poly &a, const Poly &b) {
    Poly res;
    for (auto &ka : a) {
        for (auto &kb : b) {
            TermKey tk{ka.first.b + kb.first.b, ka.first.c + kb.first.c, ka.first.d + kb.first.d};
            res[tk] += ka.second * kb.second;
        }
    }
    // remove zeros
    for (auto it=res.begin(); it!=res.end(); ) {
        if (it->second==0) it=res.erase(it); else ++it;
    }
    return res;
}

Poly poly_sub(const Poly &a, const Poly &b) {
    Poly res = a;
    poly_add(res, b, -1);
    return res;
}

Poly poly_derivative(const Poly &p) {
    Poly res;
    for (auto &kv : p) {
        long long a = kv.second;
        int b = kv.first.b, c = kv.first.c, d = kv.first.d;
        // derivative of a * x^b * sin^c * cos^d
        // term1: a*b * x^{b-1} * sin^c * cos^d
        if (b != 0) {
            TermKey tk{b-1, c, d};
            res[tk] += a * b;
        }
        // term2: a*c * x^b * sin^{c-1} * cos^{d+1}
        if (c != 0) {
            TermKey tk{b, c-1, d+1};
            res[tk] += a * c;
        }
        // term3: a*(-d) * x^b * sin^{c+1} * cos^{d-1}
        if (d != 0) {
            TermKey tk{b, c+1, d-1};
            res[tk] += -a * d;
        }
    }
    // clean zeros
    for (auto it=res.begin(); it!=res.end(); ) {
        if (it->second==0) it=res.erase(it); else ++it;
    }
    return res;
}

struct Fraction {
    Poly num, den; // denominator never zero
    Fraction() {}
    Fraction(const Poly &n, const Poly &d): num(n), den(d) {}
};

Fraction frac_add(const Fraction &f1, const Fraction &f2) {
    Poly n1 = poly_mul(f1.num, f2.den);
    Poly n2 = poly_mul(f2.num, f1.den);
    Poly num; // start with zero polynomial
    poly_add(num, n1);
    poly_add(num, n2);
    Poly den = poly_mul(f1.den, f2.den);
    return Fraction(num, den);
}

Fraction frac_sub(const Fraction &f1, const Fraction &f2) {
    Poly n1 = poly_mul(f1.num, f2.den);
    Poly n2 = poly_mul(f2.num, f1.den);
    Poly num; // start with zero polynomial
    poly_add(num, n1);
    poly_add(num, n2, -1);
    Poly den = poly_mul(f1.den, f2.den);
    return Fraction(num, den);
}

Fraction frac_mul(const Fraction &f1, const Fraction &f2) {
    Poly num = poly_mul(f1.num, f2.num);
    Poly den = poly_mul(f1.den, f2.den);
    return Fraction(num, den);
}

Fraction frac_div(const Fraction &f1, const Fraction &f2) {
    Poly num = poly_mul(f1.num, f2.den);
    Poly den = poly_mul(f1.den, f2.num);
    return Fraction(num, den);
}

Fraction frac_pow(Fraction base, long long exp) {
    Fraction res;
    // res = 1
    Poly onePoly; onePoly[{0,0,0}] = 1;
    res.num = onePoly; res.den = onePoly;
    while (exp > 0) {
        if (exp & 1) res = frac_mul(res, base);
        base = frac_mul(base, base);
        exp >>= 1;
    }
    return res;
}

Fraction frac_derivative(const Fraction &f) {
    Poly num1 = poly_mul(poly_derivative(f.num), f.den);
    Poly num2 = poly_mul(f.num, poly_derivative(f.den));
    Poly num = poly_sub(num1, num2);
    Poly den = poly_mul(f.den, f.den);
    return Fraction(num, den);
}

string term_to_string(const TermKey &tk, long long coeff, bool first) {
    string s;
    long long a = coeff;
    if (a==0) return "";
    // sign
    if (a<0) {
        s += (first?"-":"-");
        a = -a;
    } else if (!first) {
        s += "+";
    }
    // coefficient
    bool needCoeff = (a!=1) || (tk.b==0 && tk.c==0 && tk.d==0);
    if (needCoeff) s += to_string(a);
    // x part
    if (tk.b>0) {
        s += "x";
        if (tk.b>1) s += "^" + to_string(tk.b);
    }
    // sin part
    if (tk.c>0) {
        s += "sin";
        if (tk.c>1) s += "^" + to_string(tk.c);
        s += "x";
    }
    // cos part
    if (tk.d>0) {
        s += "cos";
        if (tk.d>1) s += "^" + to_string(tk.d);
        s += "x";
    }
    return s;
}

string poly_to_string(const Poly &p) {
    if (p.empty()) return "0";
    // collect keys
    vector<pair<TermKey,long long>> terms;
    for (auto &kv: p) terms.push_back(kv);
    sort(terms.begin(), terms.end(), [](auto &a, auto &b){
        if (a.first.b != b.first.b) return a.first.b > b.first.b;
        if (a.first.c != b.first.c) return a.first.c > b.first.c;
        return a.first.d > b.first.d;
    });
    string res;
    bool first=true;
    for (auto &kv: terms) {
        string part = term_to_string(kv.first, kv.second, first);
        if (!part.empty()) {
            res += part;
            first=false;
        }
    }
    if (res.empty()) return "0";
    return res;
}

string fraction_to_string(const Fraction &f) {
    // check zero numerator
    if (f.num.empty()) return "0";
    // check denominator is constant
    bool den_is_one = false;
    long long den_coeff = 0;
    if (f.den.size()==1) {
        auto it = f.den.begin();
        if (it->first.b==0 && it->first.c==0 && it->first.d==0) {
            den_coeff = it->second;
            if (den_coeff == 1) den_is_one = true;
        }
    }
    // if denominator is -1, flip numerator sign
    Poly num = f.num;
    if (!den_is_one && f.den.size()==1) {
        auto it = f.den.begin();
        if (it->first.b==0 && it->first.c==0 && it->first.d==0 && it->second == -1) {
            // multiply numerator by -1
            for (auto &kv : num) kv.second = -kv.second;
            den_is_one = true; // treat as denominator 1
        }
    }
    string numStr = poly_to_string(num);
    if (den_is_one) {
        return numStr;
    }
    string denStr = poly_to_string(f.den);
    bool num_multi = (num.size()>1);
    bool den_multi = (f.den.size()>1);
    string res;
    if (num_multi) res += "(" + numStr + ")"; else res += numStr;
    res += "/";
    if (den_multi) res += "(" + denStr + ")"; else res += denStr;
    return res;
}

// Parser
enum TokenType {TOK_END, TOK_PLUS, TOK_MINUS, TOK_MUL, TOK_DIV, TOK_POW, TOK_LP, TOK_RP, TOK_NUM, TOK_X, TOK_SIN, TOK_COS};
struct Token {TokenType type; long long value;};

struct Parser {
    string s; size_t pos; Token cur;
    Parser(const string& str): s(str), pos(0) { next(); }
    void skip_spaces(){ while(pos<s.size() && isspace(s[pos])) pos++; }
    void next(){
        skip_spaces();
        if(pos>=s.size()){ cur={TOK_END,0}; return; }
        char c=s[pos];
        if(c=='+'){ cur={TOK_PLUS,0}; pos++; }
        else if(c=='-'){ cur={TOK_MINUS,0}; pos++; }
        else if(c=='/'){ cur={TOK_DIV,0}; pos++; }
        else if(c=='*'){ cur={TOK_MUL,0}; pos++; }
        else if(c=='^'){ cur={TOK_POW,0}; pos++; }
        else if(c=='('){ cur={TOK_LP,0}; pos++; }
        else if(c==')'){ cur={TOK_RP,0}; pos++; }
        else if(isdigit(c)){
            long long val=0; while(pos<s.size() && isdigit(s[pos])){ val=val*10+(s[pos]-'0'); pos++; }
            cur={TOK_NUM,val};
        } else if(c=='x'){
            cur={TOK_X,0}; pos++; }
        else if(isalpha(c)){
            // read identifier, but handle sin/cos possibly followed by other chars (e.g., sinx)
            size_t start = pos;
            string id; while(pos<s.size() && isalpha(s[pos])){ id+=s[pos]; pos++; }
            if(id.rfind("sin",0)==0){
                cur={TOK_SIN,0};
                // consume optional following 'x'
                if(pos < s.size() && s[pos]=='x') pos++; // skip the 'x' after sin
                else pos = start + 3;
            }
            else if(id.rfind("cos",0)==0){
                cur={TOK_COS,0};
                if(pos < s.size() && s[pos]=='x') pos++; // skip the 'x' after cos
                else pos = start + 3;
            }
            else {
                // unknown identifier
                cur={TOK_END,0};
            }
        } else {
            // unknown char
            cur={TOK_END,0}; pos++; }
    }
    bool accept(TokenType t){ if(cur.type==t){ next(); return true; } return false; }
    bool expect(TokenType t){ if(cur.type==t){ next(); return true; } return false; }
    Fraction parseExpression(){
        Fraction left = parseTerm();
        while(cur.type==TOK_PLUS || cur.type==TOK_MINUS){
            bool isPlus = cur.type==TOK_PLUS; next();
            Fraction right = parseTerm();
            left = isPlus ? frac_add(left,right) : frac_sub(left,right);
        }
        return left;
    }
    // term: sequence of factors (implicit multiplication) and explicit * and / operators
    Fraction parseTerm(){
        Fraction left = parseFactor();
        while(true){
            if(cur.type==TOK_MUL){
                next();
                Fraction right = parseFactor();
                left = frac_mul(left, right);
            } else if(cur.type==TOK_DIV){
                next();
                Fraction right = parseFactor();
                left = frac_div(left, right);
            } else if(cur.type==TOK_NUM || cur.type==TOK_X || cur.type==TOK_SIN || cur.type==TOK_COS || cur.type==TOK_LP){
                // implicit multiplication
                Fraction right = parseFactor();
                left = frac_mul(left, right);
            } else break;
        }
        return left;
    }
    Fraction parseFactor(){
        // unary + or -
        if(cur.type==TOK_PLUS){ next(); return parseFactor(); }
        if(cur.type==TOK_MINUS){ next(); Fraction f = parseFactor(); // negate
            // multiply by -1 using constant -1
            Poly negNum; negNum[{0,0,0}] = -1;
            Poly one; one[{0,0,0}] = 1;
            Fraction minusOne(negNum, one);
            return frac_mul(minusOne, f);
        }
        Fraction base;
        auto make_const = [&](long long v){
            Poly num; num[{0,0,0}] = v;
            Poly den; den[{0,0,0}] = 1;
            return Fraction(num, den);
        };
        if(cur.type==TOK_NUM){ long long v=cur.value; next();
            base = make_const(v);
        } else if(cur.type==TOK_X){ next();
            Poly p; p[{1,0,0}] = 1;
            Poly one; one[{0,0,0}] = 1;
            base = Fraction(p, one);
        } else if(cur.type==TOK_SIN){
            next();
            int exp=1;
            if(cur.type==TOK_POW){ next(); if(cur.type==TOK_NUM){ exp=cur.value; next(); } }
            if(cur.type==TOK_X){ next(); }
            Poly p; p[{0,exp,0}] = 1;
            Poly one; one[{0,0,0}] = 1;
            base = Fraction(p, one);
        } else if(cur.type==TOK_COS){
            next();
            int exp=1;
            if(cur.type==TOK_POW){ next(); if(cur.type==TOK_NUM){ exp=cur.value; next(); } }
            if(cur.type==TOK_X){ next(); }
            Poly p; p[{0,0,exp}] = 1;
            Poly one; one[{0,0,0}] = 1;
            base = Fraction(p, one);
        } else if(cur.type==TOK_LP){
            next();
            base = parseExpression();
            expect(TOK_RP);
        } else {
            // error fallback to zero
            base = make_const(0);
        }
        // exponentiation
        if(cur.type==TOK_POW){
            next();
            long long exp=1;
            if(cur.type==TOK_NUM){ exp=cur.value; next(); }
            base = frac_pow(base, exp);
        }
        return base;
    }
};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    string input;
    if(!(cin>>input)) return 0;
    Parser parser(input);
    Fraction expr = parser.parseExpression();
    // normalize denominator sign if constant -1
    if (expr.den.size()==1) {
        auto it = expr.den.begin();
        if (it->first.b==0 && it->first.c==0 && it->first.d==0 && it->second==-1) {
            for (auto &kv : expr.num) kv.second = -kv.second;
            expr.den.clear(); expr.den[{0,0,0}] = 1;
        }
    }
    // output original expression as fraction (simplified) and its derivative
    string exprStr = fraction_to_string(expr);
    Fraction deriv = frac_derivative(expr);
    // normalize derivative denominator sign similarly
    if (deriv.den.size()==1) {
        auto it = deriv.den.begin();
        if (it->first.b==0 && it->first.c==0 && it->first.d==0 && it->second==-1) {
            for (auto &kv : deriv.num) kv.second = -kv.second;
            deriv.den.clear(); deriv.den[{0,0,0}] = 1;
        }
    }
    string derivStr = fraction_to_string(deriv);
    cout << exprStr << "\n" << derivStr << "\n";
    return 0;
}


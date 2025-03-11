#!/usr/bin/env python3
def curl(F):
    return f"\\nabla \\times {F}"
def derivative(F, X):
    return "\\frac{\\partial %s}{\\partial %s}" %(F, X)

def curl_vector(F, Bais = ['x', 'y', 'z']):
    n = len(Bais)
    assert(n >= 3)
    res = {
        x : "" for x in Bais
    }
    Bais = Bais * 2
    for i in range(n):
        y = Bais[i + 1]
        z = Bais[i + 2]
        res[Bais[i]] += derivative(f"{F}_{z}", y)
    for i in range(len(Bais) - 1, n - 1, -1):
        y = Bais[i - 1]
        z = Bais[i - 2]
        res[Bais[i]] += f" - {derivative(f"{F}_{z}", y)}" 
    return [v for k,v in res.items()]
def curl_to_str(vec):
    s = "("
    for i in range(len(vec)): 
        s += vec[i]
        if vec[i] != vec[-1]: s += ", "
    s +=")"
    return s
if __name__ == "__main__":
    print(curl_to_str(curl_vector("E")))
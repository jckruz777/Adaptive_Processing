import numpy as np

def aprox_inv(A, tol):
    I = np.identity(A.shape[0]) # Matriz identidad

    # Inicializacion de la matriz X
    alfa = 1 / (pow(np.linalg.norm(A, 'fro'), 2))
    X = alfa * np.conjugate(A)

    # Proceso iterativo
    while(True):
        # Condicion de parada
        error = np.linalg.norm((np.matmul(A, X) - I), 'fro')
        if (error < tol):
            break;

        X_1 = 3 * I
        X_2 = 3 * np.matmul(A, X)
        X_3 = np.matmul(A, X)
        X_3 = np.matmul(X_3, X_3)
        X = np.matmul(X, (X_1 - X_2 + X_3))

    return X

def solucion_problA(m, p, sigma, tol):
    # Calcular matriz de covarianza Rxy
    Rxx = np.ones([m, m]) * p                               # Matriz de covarianza Rxx
    np.fill_diagonal(Rxx, 1)
    A = np.eye(m,m,k=-1) + np.eye(m,m)*3 + np.eye(m,m,k=1)  # Matriz A
    At = np.conjugate(np.transpose(A))                      # Conjugado transpuesto de la matriz A
    Rxy = np.matmul(Rxx, At)

    # Calcular matriz de covarianza Ryy
    I = np.identity(m)                                      # Matriz identidad de dimension m
    Rnn = pow(sigma, 2) * I                                 # Matriz de covarianza Rnn
    Ryy = np.matmul(A, np.matmul(Rxx, At)) + Rnn            # Matriz de covarianza Ryy

    # Calular K
    K = np.matmul(Rxy, aprox_inv(Ryy, tol))                 # Matriz K del estimador de x

    # Calcular error
    error = np.trace(Rxx - np.matmul(K, np.matmul(Ryy, np.conjugate(np.transpose(K)))))

    return (K, error)

# TEST
k, error = solucion_problA(25, 0.5, 1, pow(10, -4))
print("Error: " + str(error))

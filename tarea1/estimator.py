import numpy as np
import matplotlib.pyplot as plt


def aprox_inv(A, tol):
    """
    Función que implementa la aproximación de la inversa planteado
    en el artículo 'A family of iterative methods for computing the
    approximate inverse of a square matrix and inner inverse of a
    non-square matrix'

    Parameters:
        A (Matriz): Matriz a la que se le calculará su inversa.
        tol (Tolerancia): Condición de parada del método iterativo.

    Returns:
        X: Aproximación de la inversa de la matriz A.
    """

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
    """
    Función que da solución a Problema 1 del enunciado de la
    tarea 1.

    Parameters:
        m (Dimension): Entero positivo indicando la dimension de las matrices de covarianza.
        p (rho): Valor de las entradas de la matrix de covarianza Rxx.
        sigma: Varianza que define la matriz de covarianza del ruido Rvv.
        tol (Tolerancia): Condición de parada del método iterativo.

    Returns:
        K: Matriz del estimador lineal.
        error: Error del estimador.
    """

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


def get_plot(y_min, y_max, ek_vector, x_label, y_label, title, n_points):
    """
    Función grafica el error del estimador lineal en funcion de diferentes
    variables.

    Parameters:
        y_min: Valor mínimo del eje y.
        y_max: Valor máximo del eje y.
        ek_vector: Vector con valores de error.
        x_label: Etiqueta del eje x.
        y_label: Etiqueta del eje y.
        title: Título de la gráfica.
        n_points: Número de puntos del eje y.
    """

    y = np.linspace(y_min, y_max, num=n_points)

    plt.plot(ek_vector, y, 'o-')
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()

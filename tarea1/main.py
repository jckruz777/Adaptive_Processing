import argparse
import estimator
import numpy

def main():
    # Parser de la linea de comandos
    parser = argparse.ArgumentParser(description='Un script para obtener graficas de m, sigmay rho versus el error general e_k')
    parser.add_argument('--plot', help='Tipo de grafica a obtener (m: m vs error, sigma: sigma vs error, rho: rho vs error)', default='m')
    args = parser.parse_args()
    plot_type = args.plot

    # Variables
    error_results = []
    title = ""
    x_label = ""

    # Constantes
    y_label = "Error General [e_k]"
    tol = 1e-4

    # Caso Dimension vs Error
    if plot_type == "m":
        title = "Dimension m vs Error General e_k"
        x_label = "Dimension [m]"
        print("Generando la grafica: " + title)
        rho = 0.5
        sigma = 1
        m_min = 2
        m_max = 100
        step = 1
        n_points = (m_max - m_min)/step
        for m in range(m_min, m_max, step):
            k, error = estimator.solucion_problA(m, rho, sigma, tol)
            error_results.append(error)
        estimator.get_plot(m_min, m_max, error_results, x_label, y_label, title, n_points)
    # Caso Sigma vs Error
    elif plot_type == "sigma":
        title = "Desviacion Estandar sigma vs Error General e_k"
        x_label = "Desviacion Estandar [sigma]"
        print("Generando la grafica: " + title)
        rho = 0.5
        m = 25
        sigma_min = 0.1
        sigma_max = 5.0
        step = 1
        n_points = (sigma_max - sigma_min)*10/step
        for sigma in range(int(sigma_min * 10), int(sigma_max * 10), step):
            k, error = estimator.solucion_problA(m, rho, sigma, tol)
            error_results.append(error)
        estimator.get_plot(sigma_min, sigma_max, error_results, x_label, y_label, title, n_points)
    # Caso Rho vs Error
    elif plot_type == "rho":
        title = "Constante rho vs Error General e_k"
        x_label = "Constante [rho]"
        print("Generando la grafica: " + title)
        m = 50
        sigma = 1
        rho_min = 0.01
        rho_max = 0.99
        step = 1
        n_points = (rho_max - rho_min)*100/step
        for rho in range(int(rho_min * 100), int(rho_max * 100), step):
            k, error = estimator.solucion_problA(m, rho, sigma, tol)
            error_results.append(error)
        estimator.get_plot(rho_min, rho_max, error_results, x_label, y_label, title, n_points)
    else:
        print("ERROR. Parametro invalido: " + plot_type + ".")
        print("Parametros validos: m | sigma | rho")


if __name__== "__main__":
  main()
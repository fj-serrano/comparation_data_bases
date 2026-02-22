# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import tkinter
from tkinter import filedialog

# Forma parte del import de tkinter
tk = tkinter.Tk()
tk.withdraw()  

# Función que permite al usuario seleccionar archivos
# y guardar las rutas en una lista
def select_archives():
    files = []
    name_db = []
    print("Seleccione archivos de uno en uno.\nCuando ya no quiera seleccionar mas archivos pulse cancelar.")
    while True:
        file_path = filedialog.askopenfilename()
        if file_path == "":
            break
        file_path_clean = ""
        for a in range(len(file_path)):
            if file_path[a] != "/":
                file_path_clean = file_path_clean + file_path[a]
            elif file_path[a] == "/":
                file_path_clean = file_path_clean + "\\" 
        files.append(file_path_clean)
        name_db.append(input("Nombre de la base de datos correspondiente al archivo seleccionado: "))
    return files, name_db


# Función que lee los distintos archivos seleccionados en data frames
# y extrae y opera la información importante en ellos
def files_to_df(files, name_db):
    if not files:
        print("No se han seleccionado archivos.")
        sys.exit(1)
    else:
        # Creamos la variable datos que será una lista de data frames
        datos = []
        for i in range(len(files)):
            df = pd.read_csv(files[i], sep = "\t")
            # Añadimos al data frame una nueva columna con el nombre de la base de datos
            # y otra columna con el valor absoluto del log2FoldChange para medir la desviación total
            # independientemente de si se trata de un aumento o reducción de la expresión
            df["name_db"] = name_db[i]
            df["abs_log2FoldChange"] = df["log2FoldChange"].abs()
            # Y por último en la posición de la lista datos que corresponda guardamos el data frame
            # con las columnas de interés, en este caso gene, abs_log2FoldChange y name_db
            datos.append(df[["gene", "abs_log2FoldChange", "name_db"]])
        # Una vez terminado el bucle se concatenan los data frames de la lista datos
        # en un único data frame que contendrá toda la información. Es importante activar
        # la opción de ignorar los index o nos salta un error debido a que los data frames
        # comparten index en sus columnas
        df_total = pd.concat(datos, ignore_index = True)
    return df_total
        

def main():
    files, name_db = select_archives()
    df_total = files_to_df(files, name_db)

    plt.figure()
    # Guardo en la variable db el nombre de cada base de datos
    # al almacenar todos los nombres únicos que salen
    data_base = df_total["name_db"].unique()
    contador = 0
    for db in data_base:
        # Seleccionamos la información correspondiente a la base de datos
        # que corresponde al bucle en el que estamos
        subset = df_total[df_total["name_db"] == db]
        x = np.full(len(subset), contador)
        y = subset["abs_log2FoldChange"]

        plt.scatter(x, y, alpha = 0.6)
        contador += 1

    # Ponemos los nombres a los ejes y sustituimos los valores de X por los nombres de las bases de datos
    plt.xticks(range(len(name_db)), name_db, rotation=45)
    plt.ylabel("|log2FoldChange|")
    plt.xlabel("Base de datos")
    plt.title("Comparacion de magnitud de efecto entre bases de datos")
    plt.savefig("comparation.png", bbox_inches="tight")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
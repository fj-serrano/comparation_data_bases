# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import argparse
import seaborn as sns

# FUNCIÓN QUE PERMITA AL USUARIO SELECCIONAR LOS ARCHIVOS DE ENTRADA Y GUARDAR LAS RUTAS EN UNA LISTA
def select_archives():
    files = []
    estudios = []
    # Argumentos de entrada para el script que utiliza para localizar las matrices de expresión diferencial, 
    # el nombre del estudio y localización de salida de salida
    parser = argparse.ArgumentParser(description = "Procesa un archivo de entrada")
    parser.add_argument("-i", "--input", required = True, help = "Ruta del archivo de entrada")
    parser.add_argument("-o", "--output", required = True, help = "Ruta del archivo de salida")
    parser.add_argument("-s", "--study", help = "Elige el nombre del estudio", default = "SRP")
    args = parser.parse_args()
    file_path = args.input
    output_path = args.output
    study = args.study
    # Normaliza el path para evitar problemas con rutas relativas o absolutas
    file_path_clean = os.path.normpath(os.path.abspath(file_path))
    output_path = os.path.normpath(os.path.abspath(output_path))
    # Busca todas las carpetas cuyo nombre empiece por "study" y, dentro de cada una,
    # busca recursivamente las subcarpetas de_rcadj y de_rcsa
    # Si la encuentra, añade todas las .tsv que empiecen por DESeq2, edgeR o limma a unas listas diferentes.
    metodos_prefijos = ("DESeq2", "edgeR", "limma")
    for root, dirs, files_in_dir in os.walk(file_path_clean):
        # Identifica directorios de estudio seleccionado
        if os.path.basename(root).startswith(study):
            # Si el estudio no está ya en la lista de estudios, añádelo
            if os.path.basename(root).split("_")[0] not in estudios:
                estudios.append(os.path.basename(root).split("_")[0])
            # Recorre las subcarpetas del estudio en busca de la carpeta de regresion
            for subroot, subdirs, subfiles in os.walk(root):
                # Añade solo los archivos .tsv cuyos nombres empiecen por alguno de los metodos y en los que se enfrenten 2 condiciones (vs)
                for file in subfiles:
                    if file.endswith(".tsv") and file.startswith(metodos_prefijos) and ("vs" in file):
                        file_path_full = os.path.join(subroot, file)
                        files.append(file_path_full)
    return files, output_path, estudios

# FUNCIÓN QUE LEE LOS DISTINTOS ARCHIVOS SELECCIONADOS EN DATA FRAMES, EXTRAER Y OPERAR LA INFORMACIÓN IMPORTANTE EN ELLOS
def files_to_df(files):
    if not files:
        print("No se han seleccionado archivos.")
        sys.exit(1)
    else:
        # Lista con los métodos
        metodos = ("DESeq2", "edgeR", "limma")
        # Lista con las bases de datos
        dbs = ["miRBase", "mirgenedb", "mircarta"]
        # Lee cada archivo seleccionado, cuenta el número de miRNAs significativos (padj < 0.05) y guarda esta información en un DataFrame
        # que contiene el nombre del estudio (folder_name), la base de datos usada, el método (ambos extraidos de files_names) y el número de miRNAs significativos encontrados en cada archivo
        data = []
        for file in files:
            df = pd.read_csv(file, sep = "\t")
            total_mirnas = len(df)
            significant_mirnas = df[df["padj"] < 0.05]
            significant_mirnas_log2fc = df[(df["padj"] < 0.05) & (df["log2FoldChange"].abs() > 1)]
            num_significant_mirnas = len(significant_mirnas)
            num_significant_mirnas_log2fc = len(significant_mirnas_log2fc)
            # Extrae la base de datos y el método del nombre del archivo (asumiendo que el formato es "metodo_TEXTO_DB.txt")
            metodo = "UNDEFINED"
            base_de_datos = "UNDEFINED"
            assignment = "UNDEFINED"
            # Busca el método, la base de datos, el tipo de asignacion y el nombre del estudio en la ruta del archivo
            ruta = file.split(os.sep)
            for part in ruta:
                if part.startswith("SRP"):
                    # Si empieza con SRP se asume que es el nombre del estudio, sigue la estructura "SRP_DB"
                    estudio_db = part.split("_")
                    estudio = estudio_db[0]
                    if estudio_db[1] in dbs:
                        base_de_datos = estudio_db[1]
                if part.startswith(metodos):
                    # Si empieza con método se asume que es el método, sigue la estructura "METODO_....tsv"
                    metodo_etc = part.split("_")
                    metodo = metodo_etc[0]
                if part.startswith("de_"):
                    # Si empieza con de_ se asume que es el metodo de asignación de lecturas, sigue la estructura "de_rcadj" o "de_rcsa"
                    assignment = part
            # Agrega la información a data como una lista de diccionarios
            data.append({"Estudio": estudio, "DB": base_de_datos, "Metodo": metodo, "Asignacion": assignment, "miRNAs_significativos": num_significant_mirnas, "miRNAs_significativos_DE": num_significant_mirnas_log2fc, "total_miRNAs_testeados": total_mirnas})
        # Crea un DataFrame a partir de la lista de diccionarios
        df = pd.DataFrame(data)
        df_sorted = df.sort_values(by = ["Estudio", "DB", "Metodo", "Asignacion"]).reset_index(drop = True)
    return df_sorted
      
# CREA UNA GRÁFICA DE BARRAS PARA COMPARAR EL NÚMERO DE MIRNAS SIGNIFICATIVOS ENCONTRADOS POR CADA MÉTODO Y BASE DE DATOS
def grafica_de_barras(df_total, type = "total", FC = "miRNAs_significativos", study = ""):
    # Por cada fila en el data frame, crea una barra con el número de miRNAs encontrados, 
    # se usará un color diferente para cada método
    fig, axex = plt.subplots(1, 2, figsize = (10, 8), sharey = True)
    colores = {"DESeq2": "#3182bd", "edgeR": "#e6550d", "limma": "#31a354"}

    # Primero se filtra el DataFrame para el estudio 
    df_study = df_total[df_total["Estudio"] == study].copy()
    # Y se ordena por base de datos y método
    df_study = df_study.sort_values(by = ["DB", "Metodo"])
    df_study.to_csv(f"{study}_data.csv", sep = "\t", index = False)

    # Usando los filtros se seleccionan los datos
    if type == "total":
        y_col = FC
    elif type == "porcentaje":
        df_study["porcentaje_mirnas_sig"] = (df_study[FC] / df_study["total_miRNAs_testeados"]) * 100
        y_col = "porcentaje_mirnas_sig"
    else:
        print("Tipo de grafica no reconocido")
    
    # Se obtiene la lista de asignaciones únicas para iterar sobre ellas
    asignacion = df_study["Asignacion"].unique().tolist()

    # Tras ellos iteramos y representamos cada barra con su color correspondiente
    for i, asign in enumerate(asignacion):
        ax = axex[i]
        sns.barplot(data = df_study[df_study["Asignacion"] == asign], x = "DB", y = y_col, hue = "Metodo", palette = colores, edgecolor = "white", ax = ax)
        ax.set_title(f"{asign}", fontweight = "bold", pad = 12)
        ax.legend_.remove()
        ax.set_xlabel("")
        ax.set_ylabel("")

    # Se añaden etiquetas y título a la gráfica
    title_text = type.title().replace("_", " ")
    fc_text = FC.replace("_", " ")
    fc_text = fc_text.replace("DE", "diferencialmente expresados")
    fig.suptitle(f"{title_text} de {fc_text}", fontweight = "bold")
    fig.supxlabel("Base de datos")
    fig.supylabel(f"{title_text} de {fc_text}")
    
    # Se añade una leyenda común
    handels, labels = axex[0].get_legend_handles_labels()
    fig.legend(handels[:3], labels[:3], loc = "upper left", ncol = 1, title = "Marco de Expresion Diferencial", title_fontsize = "11", fontsize = "9", bbox_to_anchor = (1.02, 1), frameon = True)

    plt.tight_layout()
    sns.despine()
    plt.savefig(f"{study}_{type}_{FC}.png", dpi = 300, bbox_inches = "tight")
    plt.close()
    return

def main():
    files, output_path, estudios = select_archives()
    df_total = files_to_df(files)
    for study in estudios:
        output_path_new = output_path + f"/{study}"
        os.makedirs(output_path_new, exist_ok = True)
        os.chdir(output_path_new)
        grafica_de_barras(df_total, study = study)
        grafica_de_barras(df_total, type = "porcentaje", study = study)
        grafica_de_barras(df_total, FC = "miRNAs_significativos_DE", study = study)
        grafica_de_barras(df_total, type = "porcentaje", FC = "miRNAs_significativos_DE", study = study)

if __name__ == "__main__":
    main()
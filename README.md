# FROM pdbmani/click_algorithm_python

Para mayores informes pueden checar el cartel de divulgacion [Detección de similitudes locales entre proteínas con dinámica molecular y teoría de gráficas](https://www.researchgate.net/publication/331327967_Deteccion_de_similitudes_locales_entre_proteinas_con_dinamica_molecular_y_teoria_de_graficas)

Se genero un script de alineamieto de dos proteinas utilizando los PDB's de las proteinas de interes, obteniendo 1 PDB alineado a la otra proteina.

Se basa en el articulo: [CLICK—topology-independent comparison of biomolecular 3D structures ](https://academic.oup.com/nar/article/39/suppl_2/W24/2506682)

Con algunas modificaciones: 

	+ filtros de angulos dihedrales
	+ filtros de ES (Estructura Secundaria)

click_align.py es el codigo que ejecuta el alineamiento le tienes que indicar la direccion de los 2 pdb's y un cutoff que consiste en el filtro de los angulos dihefrales un ejemplo seria

python Serch/click_align.py /Expermientos/pdbs/1FFT_B_1FFT_G/1FFT_B.pdb /Expermientos/pdbs/1FFT_B_1FFT_G/1FFT_G.pdb 5

Algoritmo de Reducción de Trayectorias de Plegamiento
	El algoritmo realiza la reducción de una trayectoria en dos fases. En la primera fase se forman grupos de estructuras muy cercanas tanto a nivel de posiciones de los átomos de las estructuras como a nivel de cercanía en el tiempo de aparición dentro de la trayectoria. La métrica que se utiliza en esta fase es el RMSD, que es más sensitivo a las variaciones locales de plegamiento y por lo tanto es una métrica adecuada para agrupar estructuras muy cercanas como las que aparecen de un tiempo a otro en una trayectoria de plegamiento. En la segunda fase se toman las estructuras representativas de lo grupos formados previamente y crea grupos que también tiene en cuenta la posición de los átomos entre las estructuras pero es menos sensible a las variaciones estructurales locales. La métrica para las comparaciones en esta fase es el TM-score que tiene en cuenta las propiedades globales de plegamiento y por lo tanto es adecuada para agrupar estructuras que están más alejadas o separadas tiempos más largos dentro de la trayectoria.

	El algoritmo inicia particionando la trayectoría en secciones o bins de N estructuras contiguas en el tiempo. Después sobre cada bin se realiza primero un agrupamiento local rápido y después un agrupamiento global detallado. El agrupamiento rápido aprovecha el ordenamiento temporal de las estructuras implícito en la trayectoría y forma grupos tomando inicialmente la primera estructura del bin (estado inicial) como representativa del primer grupo. La siguiente estructura del bin en orden de tiempo se compara con la representativa usando la métrica RMSD y si el valor está menor que un umbral predeterminado, la estructura se asigna a este primer grupo, de lo contrario se forma uno nuevo que toma a esta ultima estructura como su representativa. Este mismo proceso se sigue con las siguientes estructuras pero teniendo en cuenta que las comparaciones se realizan solo con las estructuras representativas de cada grupo y no con las estructuras que conforman el grupo, lo cual reduce el número de comparaciones en gran medida frente a un algoritmo convencional de agrupamiento. Además, para evitar más comparaciones, estas se realizan de atrás hacia adelante, es decir, las estructuras se comparan con la representativa del último grupo formado y si no cumple el umbral entonces se compara con el anterior a este y así sucesivamente. Esto se realiza así ya que se supone que a medida que avanza la simulación de plegamiento, las nuevas estructuras de la trayectoria van a ser más parecidas localmente (valores bajos de RMSD) a las últimas que ocurrieron. 

	En la segunda fase, el algoritmo toma las estructuras representativas de cada grupo y realiza con ellas un nuevo agrupamiento utilizando la métrica TM-score que tiene en cuenta propiedades globales del plegamiento y que no es tan sensible a variaciones locales como sucede con el RMSD. Este nuevo agrupamiento forma K grupos de los cuales se toman las estructuras centrales o medoides como representativas y que finalmente serán las que después del proceso de reducción representan a todo el bin.

	El algoritmo es fácilmente paralelizable ya que una vez particionada la trayectoría el proceso de reducción es el mismo para cada sección o bin, lo que permite que el procesamiento se reparta sobre cada bin, es decir, tanto la reducción local como la reducción global se ejecutan al mismo tiempo sobre cada bin y por lo tanto si existen N bins, cada uno de ellos se podría asignar a un proceso, hilo, o procesador.

Detalles de Implementación

	El algoritmo está implementado como a través de tres scripts: 

	• pr00_main.py: Script principal en lenguaje Python que toma los parámetros iniciales y llama a los otros scripts envíandole los parámetros necesarios.

	• pr01_createBins.py: Script en lenguaje Python que realiza la partición

	• pr02_localReduction.R : Script en lenguaje R que realiza la reducción local.

	• pr03_globalReduction.R: Script en lenguaje R que realiza la redución global..

La reducción local utiliza la función rmsd de la librería bio3d para el sistema R, y la reducción global utiliza el programa TMscore [http://cssb.biology.gatech.edu/skolnick/webservice/TM-score/index.shtml]

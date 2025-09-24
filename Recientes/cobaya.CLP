; -------------------------------------------------------------------------------------
; 1. DEFINICIÓN DE ESTRUCTURAS
; -------------------------------------------------------------------------------------

; Plantilla para almacenar los datos de la cobaya
(deftemplate cobaya
   (slot estado)
   (slot concentracion)
   (slot dosis)
   (slot condicion-fisica))

; Regla inicial para comenzar el proceso y solicitar los datos
(defrule inicio-del-programa
   (initial-fact)
   =>
   (printout t "Ingrese el estado de la enfermedad (Incipiente, Avanzado, Terminal): ")
   (bind ?estado (read))
   
   (printout t "Ingrese la concentracion de la sustancia (numero): ")
   (bind ?concentracion (read))
   
   (printout t "Ingrese el numero de dosis (numero): ")
   (bind ?dosis (read))
   
   (printout t "Ingrese la condicion fisica (Fuerte o Débil): ")
   (bind ?condicion (read))
   
   ; Afirma un nuevo hecho con los datos ingresados por el usuario
   (assert (cobaya (estado ?estado) (concentracion ?concentracion) (dosis ?dosis) (condicion-fisica ?condicion)))
)

; -------------------------------------------------------------------------------------
; 2. REGLAS DE CURACIÓN Y DEFUNCIÓN (Basadas en el análisis de la tabla)
; -------------------------------------------------------------------------------------

(defrule curacion-incipiente-debil
   (cobaya (estado "Incipiente") (condicion-fisica "Débil"))
   =>
   (printout t "Resultado: Curación." crlf))

(defrule curacion-avanzado
   (cobaya (estado "Avanzado"))
   =>
   (printout t "Resultado: Curación." crlf))

(defrule curacion-terminal-debil
   (cobaya (estado "Terminal") (condicion-fisica "Débil"))
   =>
   (printout t "Resultado: Curación." crlf))

(defrule defuncion-incipiente-concentracion
   (cobaya (estado "Incipiente") (concentracion ?c&:(> ?c 85)))
   =>
   (printout t "Resultado: Defunción." crlf))
   
(defrule defuncion-terminal-fuerte
   (cobaya (estado "Terminal") (condicion-fisica "Fuerte"))
   =>
   (printout t "Resultado: Defunción." crlf))

(defrule defuncion-por-defecto
   (cobaya (estado "Incipiente") (concentracion ?c&:(< ?c 80)))
   =>
   (printout t "Resultado: Defunción." crlf))
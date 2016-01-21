; Refraction
(define-param fcen 0.149) 
;Central fq of the Gaussian pulse
;Definition of the computational cell
(set! geometry-lattice (make lattice (size 100 50 no-size)))
; The slab
(set! geometry (list

        (make block (center -20 0) (size 1 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 2.56))))
		(make block (center -17.5 0) (size 4 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 1.63))))
		(make block (center -15 0) (size 1 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 2.56))))
		(make block (center -12.5 0) (size 4 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 1.63))))
		(make block (center -10 0) (size 1 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 2.56))))
		(make block (center -7.5 0) (size 4 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 1.63))))
		(make block (center -5 0) (size 1 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 2.56))))
		(make block (center -2.5 0) (size 4 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 1.63))))
		(make block (center 0 0) (size 1 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 2.56))))
		(make block (center 2.5 0) (size 4 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 1.63))))
		(make block (center 5 0) (size 1 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 2.56))))
		(make block (center 7.5 0) (size 4 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 1.63))))
		(make block (center 10 0) (size 1 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 2.56))))
		(make block (center 12.5 0) (size 4 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 1.63))))
		(make block (center 15 0) (size 1 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 2.56))))
		(make block (center 17.5 0) (size 4 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 1.63))))	  
		(make block (center 20 0) (size 1 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 2.56))))
		(make block (center 22.5 0) (size 4 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 1.63))))
		(make block (center 25 0) (size 1 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 2.56))))
		(make block (center 27.5 0) (size 4 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 1.63))))
		(make block (center 30 0) (size 1 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 2.56))))
		(make block (center 32.5 0) (size 4 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 1.63))))
		(make block (center 35 0) (size 1 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 2.56))))
		(make block (center 37.5 0) (size 4 infinity infinity) (e1 1 0) (e2 0 1)
                      (material (make dielectric (epsilon 1.63))))))
					  
					  
					  
					  
					  
					  
; current sources (Here,jz)
(set! sources (list
               (make source
                 (src (make continuous-src (frequency fcen)))
                 (component Ez)
                 (center -20 0)
				 (size 0 25 ))))

;PML
(set! pml-layers (list (make pml (thickness 2.0))))
                       
;Resolution
(set! resolution 5)

;Output
(use-output-directory)

; run calculation:
(run-until 300  
               (at-beginning output-epsilon)
               (to-appended "ez" (at-every 0.6 output-efield-z)))

           
	       
       
	
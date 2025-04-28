#lang racket

(require plot)
(require math/matrix)
(require racket/random)
(random-seed 69)
(provide (all-defined-out))


#|

  Nonlinearity with Basis Expansions

  All below code is trying to regress data using transformation in for of:

  f(X) - ∑β_m·h_m(X)     with sum over m = 1,2,..., M

  In general each h_m(X) is transformation ℝ^p -> ℝ, but in this code below I'm only doing p=1.

|#


(define-values (DATA_MIN_X ξ1 ξ2 DATA_MAX_X) (values 0.0 1.0 2.0 3.0))


#| DATA |#
(define (data-fun x)
  (/ (sin (* 3 x))
     4))

(define data
  (for/list ([x (range DATA_MIN_X DATA_MAX_X 0.1)])
    (list x (data-fun x))))



;; generate one N(0,1) sample
(define (gaussian-standard)
  (let ([u1 (random)]
        [u2 (random)])
    (* (sqrt (* -2 (log u1)))
       (cos (* 2 pi u2)))))

;; give it mean μ and sd σ
(define (gaussian-noise mean sigma)
  (+ mean (* sigma (gaussian-standard))))

(define noisy-data
  (map (λ (pt)
         (list (first pt)
               (gaussian-noise (second pt) 0.1)))
       data))




#| Piecewise Constant |#
(define (h-interval #:left-limit [left-limit #f]
                    #:right-limit [right-limit #f])
  (lambda (x)
    (if (and (if right-limit
                 (< x right-limit)
                 #t)
             (if left-limit
                 (>= x left-limit)
                 #t))
        1
        0)))

(define H-constant
  (list (h-interval #:right-limit ξ1)
        (h-interval #:left-limit ξ1 #:right-limit ξ2)
        (h-interval #:left-limit ξ2)))

(define (β-constant X
            #:left-limit [left-limit #f]
            #:right-limit [right-limit #f])
  (let ([out (filter (λ (x) (not (false? x)))
                     (for/list ([x X])
                       (if (and (if right-limit ; no right limit if not defined
                                    (< (first x) right-limit)
                                    #t)
                                (if left-limit  ; not left limit if not defined
                                    (>= (first x) left-limit)
                                    #t))
                           (second x)
                           #f)))])
    (if (> (length out) 0)
        (/ (foldl + 0 out) (length out))
        0)))


(define β-piecewise-constant
  (list (β-constant noisy-data #:right-limit ξ1)
        (β-constant noisy-data #:left-limit ξ1 #:right-limit ξ2)
        (β-constant noisy-data #:left-limit ξ2)))



(define (piecewise-constant x)
  (foldl +
         0
         (map (λ (βi hi)
                (* βi
                   (hi x)))
              β-piecewise-constant
              H-constant)))





#| Piecewise Linear |#
(define H-linear
  (list (h-interval #:right-limit ξ1)
        (h-interval #:left-limit ξ1 #:right-limit ξ2)
        (h-interval #:left-limit ξ2)))



(define (β-piecewise-linear #:left-limit [left-limit DATA_MIN_X]
                            #:right-limit [right-limit DATA_MAX_X])
  (let ([relevant-data-range (filter (λ (x) (not (false? x)))
                                     (for/list ([x noisy-data])
                                       (if (and (if right-limit ; no left limit if not defined
                                                    (< (first x) right-limit)
                                                    #t)
                                                (if left-limit  ; not right limit if not defined
                                                    (>= (first x) left-limit)
                                                    #t))
                                           x
                                           #f)))])
    (let ([X (list*->matrix (map (λ (data) (list 1 (first data)))
                                                               relevant-data-range))]
          [Y (list*->matrix (map (λ (data) (cdr data))
                                                               relevant-data-range))])
      (matrix->list
       (matrix* (matrix-inverse (matrix* (matrix-transpose X) X))
                (matrix* (matrix-transpose X) Y))))))


(define β-linear
  (list (β-piecewise-linear #:right-limit ξ1)
        (β-piecewise-linear #:left-limit ξ1 #:right-limit ξ2)
        (β-piecewise-linear #:left-limit ξ2)))


(define (piecewise-linear x)
  (foldr +
         0
         (map (λ (βi hi)
                (+
                 (* (car βi)
                    (hi x))
                 (* (car (cdr βi))
                   (* x
                      (hi x)))))
              (append β-linear)
              (append H-linear))))











; Continuous piecewise linear with hand-picked parameters, NOT OPTIMAL
(define H-continuous
  (list (λ (x) 1)
        (λ (x) x)
        (λ (x) (let ([out (- x ξ1)])
                 (if (> out 0)
                     out
                     0)))
        (λ (x) (let ([out (- x ξ2)])
                 (if (> out 0)
                     out
                     0)))))

(define (continuous-piecewise-linear x)
  (foldr +
         0
         (map (λ (βi hi)
                (* βi
                   (hi x)))
              (list 0.3 -0.3 -0.05 1.4)
              H-continuous)))









#| Cubic Spline

  Splines are special case of basis function regressions.
  Their basis funcitons correspond to polynomial terms and has additional terms for each knot.
  For instance Spline of order M=4 with K=2 knots will have M+K=6 terms:
  1. "Monomial terms" (it's not official name, I just call them that way for understanding):
    h_1(x) = 1
    h_2(x) = x
    h_3(x) = x^2
    h_4(x) = x^3
  2. "Knot terms" (naming as above):
    h_5(x) = RELU((x-ξ_1)^3)
    h_6(x) = RELU((x-ξ_2)^3)

  The "Knot" terms are always degree M, so you can see they are basically the same as h_M,
                                        but moved right by ξ and with additional RELU on top.

  Spline of degree 4 (M=4) is commonly named Cubic Spline.
|#
(define (piecewise-spline-4 β)
  (define H
    (list (λ (x) 1)
          (λ (x) x)
          (λ (x) (expt x 2))
          (λ (x) (expt x 3))
          (λ (x) (let ([out (expt (- x ξ1) 3)])
                   (if (> out 0)
                       out
                       0)))
          (λ (x) (let ([out (expt (- x ξ2) 3)])
                   (if (> out 0)
                       out
                       0)))))

  (λ (x) (foldr +
                0
                (map (λ (βi hi)
                       (* βi
                          (hi x)))
                     β
                     H))))

#|
  The function above defines Spline function, but β coefficients must be defined.
  There is a closed-form solution:
  β^ = (HX^T*HX + λ*Ω)^{-1} * (HX^T * Y)

  where Ω[j,k] = ∫ h_j''(x) * h_k''(x) dx     , j = 1, 2, ..., M+K
                                                k = 1, 2, ..., M+K
        HX[i,j] = h_j(x_i)                    , i = 1, 2, ..., N
                                                j = 1, 2, ..., M+K
        Y[i] = y_i                            , i = 1, 2, ..., N

  And λ controls smoothness of the curve, it's penalty component known from Ridge Regression,
    but here penalty is different.
|#



; matrix Nx(M+K) in form of: HX_{i,j} = h_j(x_i)
;         where h_j are spline basis functions j∈0, 1, ..., M+K-1
;         and M=4 (cubic spline) and K=2 (number of knots: |{ξ1, ξ2}| = 2)
(define HX
  [build-matrix (length noisy-data)
                6
                (λ (i j) (let ([x (car (list-ref noisy-data i))])
                           (match j
                           [0 1]
                           [1 x]
                           [2 (expt x 2)]
                           [3 (expt x 3)]
                           [4 (let ([out (expt (- x ξ1) 3)])
                                (if (> out 0)
                                    out
                                    0))]
                           [5 (let ([out (expt (- x ξ2) 3)])
                                (if (> out 0)
                                    out
                                    0))]
                           [else (error "wrong id")])))])

; real/reference data output
(define Y
  [build-matrix (length noisy-data)
                1
                (λ (i j) (second (list-ref noisy-data i)))])

;; Omega[j,k] = ∫ h_j''(x) * h_k''(x) dx
;; h'' is: [0 0 2 6x 6(x-ξ1) 6(x-ξ2)]
#;(define (Ω-antiderivative x)
  (matrix
   [[0 0 0          0                    0                                       0                                     ]
    [0 0 0          0                    0                                       0                                     ]
    [0 0 4x         12x                  12(x-ξ1)                                12(x-ξ2)                              ]
    [0 0 12x        12x^3                (12x^3 - 18ξ1*x^2)                      (12x^3 - 18ξ2*x^2)                    ]
    [0 0 12(x-ξ1)   (12x^3 - 18ξ1*x^2)   (12x^3 - 36ξ1*x^2 + 36ξ1^2*x)           (12x^3 - 18(ξ1 + ξ2)*x^2 + 36ξ1*ξ2*x) ]
    [0 0 12(x-ξ2)   (12x^3 - 18ξ2*x^2)   (12x^3 - 18(ξ1 + ξ2)*x^2 + 36ξ1*ξ2*x)   (12x^3 - 36ξ2*x^2 + 36ξ2^2*x)         ]]))
;; Same as above but in prefix notation:
(define (Ω-antiderivative x)
  (matrix
   [[0 0 0 0 0 0]
    [0 0 0 0 0 0]
    [0 0 (* 4 x)                      (* 12 x)                   (* 12 (- x ξ1))                  (* 12 (- x ξ2))]
    [0 0 (* 12 x)                     (* 12 (expt x 3))          (- (* 12 (expt x 3))
                                                                    (* 18 ξ1 (expt x 2)))         (- (* 12 (expt x 3))
                                                                                                     (* 18 ξ2 (expt x 2)))]
    [0 0 (* 12 (- x ξ1))              (- (* 12 (expt x 3))
                                         (* 18 ξ1 (expt x 2)))   (+ (- (* 12
                                                                          (expt x 3))
                                                                       (* 36 ξ1
                                                                          (expt x 2)))
                                                                    (* 36 (expt ξ1 2) x))         (+ (- (* 12 (expt x 3))
                                                                                                        (* 18 (+ ξ1 ξ2)
                                                                                                           (expt x 2)))
                                                                                                     (* 36 ξ1 ξ2 x))]
   [0 0 (* 12 (- x ξ2))               (- (* 12 (expt x 3))
                                         (* 18 ξ2 (expt x 2)))   (+ (- (* 12 (expt x 3))
                                                                       (* 18 (+ ξ1 ξ2)
                                                                          (expt x 2)))
                                                                    (* 36 ξ1 ξ2 x))               (+ (- (* 12 (expt x 3))
                                                                                                        (* 36 ξ2 (expt x 2)))
                                                                                                     (* 36 (expt ξ2 2) x))]]))

(define Ω
  (matrix- (Ω-antiderivative DATA_MAX_X)
           (Ω-antiderivative DATA_MIN_X)))


(define (β^ lambda_)
  (matrix->list
   (matrix* (matrix-inverse
             (matrix+ (matrix* (matrix-transpose HX)
                               HX)
                      (matrix-scale Ω lambda_)))
            (matrix* (matrix-transpose HX)
                     Y))))


(define spline-4 (piecewise-spline-4 (β^ 0.0)))

(cons
 (plot
  (list
   (function data-fun #:color "blue" #:label "True Function")
   (points noisy-data #:label "Real Data")
   (function piecewise-constant #:color "green" #:label "Piecewise Constant")
   (function piecewise-linear #:color "violet" #:label "Piecewise Linear")
   (function continuous-piecewise-linear #:color "orange" #:label "Continuous Piecewise Linear")
   (vrule ξ1 #:style 'long-dash #:label "ξ1")
   (vrule ξ2 #:style 'long-dash #:label "ξ2"))
  #:y-min -0.7
  #:y-max 0.7)

 (plot
 (list
  (function data-fun #:color "blue" #:label "True Function")
  (points noisy-data #:label "Real Data")
  (function spline-4 #:color "red" #:label "Spline (M=4 aka Cubic)")
  (vrule ξ1 #:style 'long-dash #:label "ξ1")
  (vrule ξ2 #:style 'long-dash #:label "ξ2"))
 #:y-min -0.7
 #:y-max 0.7))






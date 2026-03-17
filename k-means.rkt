#lang racket

(require plot)


(define (k-means datapoints number-centroids
                 #:initial-centroids-positions [initial-centroids-positions #f]
                 #:min-change-required [min-change-required 0.001]
                 #:max-iterations [max-iterations 50])
  (begin
    (define num-dimensions (vector-length (car datapoints)))

    
    (define centroids
      (if (not initial-centroids-positions)
          (for/list ([i (range number-centroids)])
            (build-vector num-dimensions (λ (_) (random))))
          (if (not (equal? (length initial-centroids-positions)
                          number-centroids))
              (error
               (format "length of initial-centroids-positions (~a) must be equal to number-centroids (~a), but wasn't"
                       (length initial-centroids-positions)
                       number-centroids))
              initial-centroids-positions)))

    
    (define (datapoints-with-centroids centroids)
      (foldl
       (λ (datapoint datapoints-acc)
         (cons 
               (cons datapoint
                     (car
                      (argmin cdr
                              (for/list ([(centroid i) (in-indexed centroids)])
                                (cons i
                                      (foldl
                                       (λ (dim acc)
                                         (+ acc
                                            (expt (- (vector-ref centroid dim)
                                                     (vector-ref datapoint dim))
                                                  2)))
                                       0
                                       (range num-dimensions)))))))
               datapoints-acc))
       '()
       datapoints))

    
    (define previous-centroids #f)
        
    (for ([i (in-range max-iterations)]
                #:unless (and previous-centroids
                              (< (for/fold ([total 0])
                                           ([previous previous-centroids]
                                            [this centroids])
                                   (sqrt
                                    (foldl
                                     (λ (dim acc)
                                       (+ acc
                                          (expt (- (vector-ref previous dim)
                                                   (vector-ref this dim))
                                                2)))
                                     0 (range num-dimensions))))
                                 min-change-required))
                )
      (begin
        (set! previous-centroids centroids)
        (set! centroids
              (let ([dp-with-centroids (datapoints-with-centroids centroids)])
                (for/foldr ([acc '()])
                  ([centroid-id (in-range (length centroids))])
                  (let ([datapoints-in-centroid (map car (filter (λ (datapoint-with-centroid)
                                                                   (equal? (cdr datapoint-with-centroid)
                                                                           centroid-id))
                                                                 dp-with-centroids))])
                    (cons (if (> (length datapoints-in-centroid) 0)

                              (vector-map (λ (x) (/ x (length datapoints-in-centroid)))
                                          (apply vector-map + 
                                                 datapoints-in-centroid))
                              
                              (list-ref centroids centroid-id))
                          acc)))))))

    (values centroids (datapoints-with-centroids centroids))))


(define datapoints
  (list (vector 0.2 0.3)
        (vector 0.5 0.5)
        (vector 0.7 0.1)))


#;(define datapoints
  (list (vector 0.2)
        (vector 0.5)
        (vector 0.7)))

(define-values (centroids datapoints-with-centroids)
  (k-means datapoints 2))


(define (pick-centroid . dim-components) ;; x, y, z...
  (cdr
   (argmin car (map (λ (centroid)
                      (cons
                       (for/fold ([acc 0])
                                 ([dim (min (length dim-components)
                                            (vector-length centroid))])
                         
                         (+ acc
                            (expt (- (vector-ref centroid dim)
                                     (list-ref dim-components dim))
                                  2)))
                            centroid))
                    centroids))))


(define (k-means-vector x y)
  (let ([centroid-point (pick-centroid x y)])
    (match (vector-length centroid-point)
      [0 (error "centroid-point has shape of 0... not acceptible")]
      [1 (vector (- (vector-ref centroid-point 0) x)
                 0)]
      [_ (vector (- (vector-ref centroid-point 0) x)
                 (- (vector-ref centroid-point 1) y))])))



(define (plot-datapoints datapoints centroids)
  (let ([dims (vector-length (car datapoints))])
    (match dims
      [0 (error "centroid-point has shape of 0... not acceptible")]
      [1 (plot
           (list (vector-field k-means-vector
                         -0.1 1.1 -0.5 0.5
                         #:alpha 0.4)
                
                 (for/list ([datapoint datapoints-with-centroids])
                   (begin
                     (define colors '("red" "darkblue" "darkgreen" "orange" "purple" "cyan" "magenta"
                                            "brown" "pink" "lime" "teal" "navy" "maroon" "olive" "coral"))
                     (points (map (λ (dp) (vector (vector-ref dp 0) 0)) (list (car datapoint)))
                             #:color "darkblue"
                             #:fill-color (list-ref colors (cdr datapoint))
                             #:sym 'fullcircle5
                             #:size 5)))
                 #;(points (map (λ (c) (vector (vector-ref c 0) 0)) centroids)
                           #:color "darkgreen"
                           #:size 30))
          #:height 300
          #:width 1500
          #:x-min -0.1
          #:x-max 1.1
          #:y-min -0.5
          #:y-max 0.5)]
      [_ (plot

          (list (vector-field k-means-vector
                              0 1 0 1
                              #:alpha 0.4)
                
                (for/list ([datapoint datapoints-with-centroids])
                  (begin
                    (define colors '("red" "darkblue" "darkgreen" "orange" "purple" "cyan" "magenta"
                                           "brown" "pink" "lime" "teal" "navy" "maroon" "olive" "coral"))
                    (points (map (λ (dp) (vector-copy dp 0 2)) (list (car datapoint)))
                          #:color "darkblue"
                          #:fill-color (list-ref colors (cdr datapoint))
                          #:sym 'fullcircle5
                          #:size 10)))
                
                #;(points (map (λ (c) (vector-copy c 0 2)) centroids)
                   #:color "darkgreen"
                   #:size 30))
          #:height 800
          #:width 800
          #:x-min -0.1
          #:x-max 1.1
          #:y-min -0.1
          #:y-max 1.1)])))


(plot-datapoints datapoints centroids)





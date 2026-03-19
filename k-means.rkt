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
                #:break (and previous-centroids
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


(define (generate datapoints)

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


  ;(println datapoints-with-centroids)

  (define (k-means-vector x y)
    (let ([centroid-point (pick-centroid x y)])
      (match (vector-length centroid-point)
        [0 (error "centroid-point has shape of 0... not acceptible")]
        [1 (vector (- (vector-ref centroid-point 0) x)
                   0)]
        [_ (vector (- (vector-ref centroid-point 0) x)
                   (- (vector-ref centroid-point 1) y))])))





  (define (sillhouette-scores datapoints-with-centroids)
    (begin
      (define closest-cluster-distance
        (for/fold ([acc (make-hash)])
                  ([datapoint datapoints-with-centroids])
          (hash-set! acc (car datapoint)
                     (make-hash
                      (list (cons 'cluster-id (cdr datapoint))
                            (cons 'in-group-cluster
                                  (make-hash '((distances . ()) (score . 0))))
                            (cons 'out-group-clusters
                                  (for/fold ([acc-clusters (make-hash)])
                                            ([cluster-id (in-range (length centroids))]
                                             #:unless (equal? (cdr datapoint) cluster-id))
                                    (hash-set! acc-clusters cluster-id
                                               (make-hash (list (cons 'distances '())
                                                                (cons 'score 0))))
                                    acc-clusters)))))
          acc))

      (for/list ([(datapoint-left l-id) (in-indexed datapoints-with-centroids)])
        (for/list ([(datapoint-right r-id) (in-indexed datapoints-with-centroids)]
                   #:unless (equal? l-id r-id))
          (let
              ([distance (sqrt (apply +
                                      (vector->list
                                       (vector-map (λ (vec) (expt vec 2))
                                                   (vector-map (λ (vec1 vec2) (- vec1 vec2))
                                                               (car datapoint-left)
                                                               (car datapoint-right))))))])
               
            (if (equal? (cdr datapoint-left) (cdr datapoint-right))
                   
                (let ([new-distances
                       (append
                        (hash-ref (hash-ref (hash-ref closest-cluster-distance
                                                      (car datapoint-left))
                                            'in-group-cluster)
                                  'distances)
                        (list distance))])
                  (hash-set! (hash-ref (hash-ref closest-cluster-distance
                                                 (car datapoint-left))
                                       'in-group-cluster)
                             'distances
                             new-distances))
                   
                (let ([new-distances
                       (append
                        (hash-ref (hash-ref (hash-ref (hash-ref closest-cluster-distance
                                                                (car datapoint-left))
                                                      'out-group-clusters)
                                            (cdr datapoint-right))
                                  'distances)
                        (list distance))])
                  (hash-set! (hash-ref (hash-ref (hash-ref closest-cluster-distance
                                                           (car datapoint-left))
                                                 'out-group-clusters)
                                       (cdr datapoint-right))
                             'distances
                             new-distances)))
            distance)))

      
      (for/list ([(datapoint data) closest-cluster-distance])
        (let ([a (let ([len (length (hash-ref (hash-ref data 'in-group-cluster) 'distances))])
                   (match len
                     [0 0]
                     [1 (car (hash-ref (hash-ref data 'in-group-cluster) 'distances))]
                     [_ (/ (apply + (hash-ref (hash-ref data 'in-group-cluster) 'distances))
                           len)]))]
              [b (foldl min 99999
                        (for/list ([(cluster-id data) (hash-ref data 'out-group-clusters)]
                                   #:unless (equal? 0 (length (hash-ref data 'distances))))
                          (let ([len (length (hash-ref data 'distances))])
                            (match len
                              [0 0]
                              [1 (car (hash-ref data 'distances))]
                              [_ (/ (apply + (hash-ref data 'distances))
                                    len)]))))])
          (list datapoint
                (hash-ref data 'cluster-id)
                (/ (- b a)
                   (max a b)))))))

      
  (values (sillhouette-scores datapoints-with-centroids) centroids k-means-vector))




(define (plot-datapoints sillhuette-scores centroids k-means-vector
                         #:with-vec-field [with-vec-field #f]
                         #:with-centroid-points [with-centroid-points #t]
                         #:with-arrows [with-arrows #t])
  (define colors '("red" "darkblue" "darkgreen" "orange" "purple" "cyan" "magenta"
                         "brown" "pink" "lime" "teal" "navy" "maroon" "olive" "coral"))
  (let ([dims (vector-length (car (car sillhuette-scores)))])
    (match dims
      [0 (error "centroid-point has shape of 0... not acceptible")]
      [1 (plot
          (list (if with-vec-field
                    (vector-field k-means-vector
                                  -0.1 1.1 -0.5 0.5
                                  #:alpha 0.4)
                    '())
                 
                (for/list ([datapoint sillhuette-scores])
                  (point-label (vector (vector-ref (car datapoint) 0) 0)
                               (real->decimal-string (car (cddr datapoint)) 2)
                               #:point-color (list-ref colors (modulo (car (cdr datapoint)) (length colors)))
                               #:point-fill-color (list-ref colors (modulo (car (cdr datapoint)) (length colors)))
                               #:point-sym 'fullcircle5
                               #:color (list-ref colors (modulo (car (cdr datapoint)) (length colors)))
                               #:point-size 10
                               #:size 20))

                (if with-arrows
                    (for/list ([datapoint sillhuette-scores]
                               #:when (> (abs (- (vector-ref (car datapoint) 0)
                                                 (vector-ref (list-ref centroids (car (cdr datapoint))) 0)))
                                         0))
                      (arrows (list (vector (vector-ref (car datapoint) 0) 0)
                                    (vector (vector-ref (list-ref centroids (car (cdr datapoint))) 0) 0))
                              #:arrow-head-size-or-scale '(= 10)))
                    '())
                 
                (if with-centroid-points
                    (points (map (λ (c) (vector (vector-ref c 0) 0)) centroids)
                            #:color "darkgreen"
                            #:sym 'fullcircle5
                            #:size 3)
                    '()))
          #:height 300
          #:width 1500
          #:x-min -0.1
          #:x-max 1.1
          #:y-min -0.5
          #:y-max 0.5)]
      [_ (plot

          (list (if with-vec-field
                    (vector-field k-means-vector
                                  0 1 0 1
                                  #:alpha 0.4)
                    '())

                (for/list ([datapoint sillhuette-scores])
                  (point-label (vector-copy (car datapoint) 0 2)
                               (real->decimal-string (car (cddr datapoint)) 2)
                               #:point-color (list-ref colors (modulo (car (cdr datapoint)) (length colors)))
                               #:point-fill-color (list-ref colors (modulo (car (cdr datapoint)) (length colors)))
                               #:point-sym 'fullcircle5
                               #:color (list-ref colors (modulo (car (cdr datapoint)) (length colors)))
                               #:point-size 10
                               #:size 20))

                (if with-arrows
                    (for/list ([datapoint sillhuette-scores]
                               #:when (> (abs (- (vector-ref (car datapoint) 0)
                                                 (vector-ref (list-ref centroids (car (cdr datapoint))) 0)))
                                         0))
                      (arrows (list (vector-copy (car datapoint) 0 2)
                                    (vector-copy (list-ref centroids (car (cdr datapoint))) 0 2))
                              #:arrow-head-size-or-scale '(= 10)))
                    '())

                (if with-centroid-points
                    (points (map (λ (c) (vector-copy c 0 2)) centroids)
                            #:color "darkgreen"
                            #:sym 'fullcircle5
                            #:size 3)
                    '()))
          #:height 800
          #:width 800
          #:x-min -0.1
          #:x-max 1.1
          #:y-min -0.1
          #:y-max 1.1)])))



(define datapoints
  (list (vector 0.2 0.3)
        (vector 0.5 0.5)
        (vector 0.7 0.1)))


#;(define datapoints
    (list (vector 0.2)
          (vector 0.5)
          (vector 0.7)))

(define-values (sillhuette-scores centroids k-means-vector) (generate datapoints))
;(println sillhuette-scores)
(plot-datapoints sillhuette-scores centroids k-means-vector)





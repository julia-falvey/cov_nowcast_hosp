
pull_api = function(endpoint, limit = "5000000") {
  return(read_csv(paste0(endpoint, "?$limit=", limit)))
}

scrape_list = c('https://data.cdc.gov/resource/g653-rqe2.csv', 
                'https://data.cdc.gov/resource/2ew6-ywp6.csv',
                'https://data.cdc.gov/resource/rdmq-nq56.csv',
                'https://data.cdc.gov/resource/seuz-s2cv.csv', 
                'https://data.cdc.gov/resource/d2tw-32xv.csv', 
                'https://data.cdc.gov/resource/275g-9x8h.csv',
                'https://data.cdc.gov/resource/i2a4-xk9k.csv',
                'https://data.cdc.gov/resource/kipu-qxy8.csv',
                'https://data.cdc.gov/resource/twtx-bfcw.csv', 
                'https://data.cdc.gov/resource/cf5u-bm9w.csv',
                'https://data.cdc.gov/resource/bigw-pgk2.csv',
                'https://data.cdc.gov/resource/gvsb-yw6g.csv',
                'https://data.cdc.gov/resource/jr58-6ysp.csv',
                'https://data.cdc.gov/resource/9t9r-e5a3.csv',
                'https://data.cdc.gov/resource/7xva-uux8.csv',
                'https://data.cdc.gov/resource/n8mc-b4w4.csv',
                'https://data.cdc.gov/resource/n8mc-b4w4.csv',
                'https://data.cdc.gov/resource/kvib-3txy.csv',
                'https://data.cdc.gov/resource/7dk4-g6vg.csv',
                'https://data.cdc.gov/resource/39z2-9zu6.csv')

full_scrape = lapply(scrape_list, function(x) pull_api(x))
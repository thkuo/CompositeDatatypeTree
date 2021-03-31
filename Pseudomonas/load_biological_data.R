## import the clinical data of the samples
load_biological_data<-function(){
    bio_df_f<- './hospitals.csv'
    bio_df<- read.csv(bio_df_f, header= T, quote = '"' , check.names = F,
                      stringsAsFactors = F)
    
    col_names<- colnames(bio_df)
    
    bio_df$Isolate<- gsub(' ', '', bio_df$Isolate)
    rownames(bio_df)<- bio_df$Isolate
    colnames(bio_df)<- col_names   
    bio_df
}

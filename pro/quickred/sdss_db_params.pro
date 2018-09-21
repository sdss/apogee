pro SDSS_DB_PARAMS

;----------------------------------------------------
;
;  IDL startup file for the IDL-SQL tool to
;  connect to the SDSS APO database
;
;----------------------------------------------------
;
    ; Information for the SDSS4 production database
    defsysv, '!_IDL_SQL_DRIVER','org.postgresql.Driver' 
    defsysv, '!_IDL_SQL_USER','sdssdb_admin'
    defsysv, '!_IDL_SQL_PASS','4-photometry'
    defsysv, '!_IDL_SQL_PROTOCOL','jdbc:postgresql://'
    defsysv, '!_IDL_SQL_HOST','sdss4-db'
    defsysv, '!_IDL_SQL_PORT','5432'
    defsysv, '!_IDL_SQL_DB','apodb'

    ;defsysv, '!_IDL_SQL_DRIVER','org.postgresql.Driver' 
    ;defsysv, '!_IDL_SQL_USER','postgres'
    ;defsysv, '!_IDL_SQL_PASS',''
    ;defsysv, '!_IDL_SQL_PROTOCOL','jdbc:postgresql://'
    ;defsysv, '!_IDL_SQL_HOST','sdss-db-p'
    ;defsysv, '!_IDL_SQL_PORT','5432'
    ;defsysv, '!_IDL_SQL_DB','apo_platedb'
end

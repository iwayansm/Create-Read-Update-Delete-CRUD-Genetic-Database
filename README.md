# Create-Read-Update-Delete-CRUD-Genetic-Database
PROJECT TITLE:
THE NATIONAL GENETIC DATABASE SYSTEM OF INDONESIA
------THIS IS A PART OF CAPSTONE PROJECT MODULE-1 PURWADHIKA------

The objective is to construct a program capable of performing the four fundamental operations of data management: creation, reading, updating, and deletion (CRUD). These operations pertain to the genetic profile of species of organisms, inclusive of the DNA sequence. The program utilise the regular functions within the Python programming language.
In this project, I have elected to develop a program that enables the management of the genetic record in the database through the execution of six distinct operations. These are: four CRUD (Create, Read, Update, Delete) operations; a singular user registration and sign-in operation; and an operation to initiate, modify, and save the database.

The functions are as follows:

A.  SINGULAR USER REGISTRATION AND SIGN-IN
    The system permits the user to register their identity and 16-digit identity card number, as in the case of the 'Kartu Tanda Penduduk' (KTP). Upon completion of the registration process, the user is then able
    to sign in using the registered research ID and ID card number.
    DISCLAIMER: At the time of writing, the system is not yet connected to any known servers that store the social ID card data of Indonesia. At this time, neither the identity number nor the identity card number
                are stored in any file or other repository.

B.  CREATE-READ-UPDATE-DELETE (CRUD)
    1.  SUBMIT RECORD (CREATE)
        The system enables the user to submit multiple, non-duplicated records of genetic profiles, which include the species name. At present, the system permits the submission of DNA sequences of a maximum
        length of 33 base pairs. Upon completion of the submission process, a submission metadata of the record will be generated, along with unique identifier or Reference ID (ID.Ref), which must be registered
        in subsequent operations.
    2.  FETCH RECORD (READ)
        The system permits registered users to view the entirety of registered records in the genetic database, or a selected subset thereof, by inputting the unique identifier (ID.Ref) generated during the
        submission process.
    3.  REVIEW RECORD (UPDATE)
        The system enables the user to review the record in the database by inputting the unique identifier (ID.Ref) from the submission operation. The user is permitted to modify the data pertaining to
        the species name, date, location, sequencing methodology, and DNA sequence. Upon altering the species name, date, location, or sequencing methodology, a new metadata record and a unique identifier (ID.Ref)
        will be generated, replacing the existing metadata record.
    4.  WITHDRAW RECORD (DELETE)
        The system enables the user to retrieve a particular record from the genetic database, based on the ID.Ref from the metadata generated during the submission process.

C.  INITIATE-MODIFY-SAVE THE DATABASE
    The genetic database is initiated by utilising a pickle file as the handle for the execution of the CRUD operation. The information pertaining to the record and the entirety of the genetic database is
    essentially stored in a pickle file, which enables the user to create a new file, store information in the file, modify information in the file, and save the information in the file. Prior to the commencement
    of the process, the record will be in its original, unpopulated state. The function will then generate a new, empty file. Once the CRUD operation has been completed, any alterations will be saved prior to
    exiting the program.
     


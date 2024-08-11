# Create-Read-Update-Delete-CRUD-Genetic-Database
PROJECT THEME: DATA SCIENCE

PROJECT TITLE: THE NATIONAL GENETIC DATABASE SYSTEM OF INDONESIA

------THIS IS A PART OF CAPSTONE PROJECT MODULE-1 PURWADHIKA------

The objective is to create a program that can perform the four fundamental data management operations: create, read, update, and delete (CRUD). These operations are related to the collections of DNA sequences of various species. The program uses regular functions within the Python programming language environment in Visual Studio Code.

In this project, I am trying to develop a prototype program to manage genetic records in a database. The program will perform six main operations: four CRUD (Create, Read, Update, Delete) operations, a user registration and sign-in operation, and a database initialisation, modification, and saving operations.

The functions are listed below:

A. USER REGISTRATION AND SIGN-IN

The system enables users to register their identity with their 16-digit identity card number as a password, as listed in the Indonesian social ID card, known locally as 'Kartu Tanda Penduduk' or KTP. After completing the registration process, the user can sign in using the registered research ID and ID card number.

**DISCLAIMERS**
The system is currently not connected to any known servers storing social ID card data of Indonesia, nor is it stored in any read-only file or other repository.

B. CREATE-READ-UPDATE-DELETE (CRUD)

1.  SUBMIT RECORD (CREATE)

    The system allows users to submit multiple non-duplicated records
	of genetic profiles that consist of species name, date, location,
	sequencing methodology, and DNA sequence. Currently, the system
	only accepts DNA sequences up to a maximum length of 33 base
	pairs. After the submission process is complete, metadata for the
	submission will be generated, along with a unique identifier or
	Reference ID (ID.Ref), which is required to access other functions
	and further operations.

2.  FETCH RECORD (READ)

    The system allows registered users to access all registered records
	in the genetic database or a selected subset by entering the ID.Ref
        obtained during the submission process.

3.  REVIEW RECORD (UPDATE)

    To review a record in the database, the user can input the ID.Ref
	from the submission operation. Then, the user can modify the 
	species name, date, location, sequencing methodology, and DNA
	sequence. If the user changes any of these details, a new
	metadata record and an ID.Ref will be generated to replace the
	existing record.

4.  WITHDRAW RECORD (DELETE)
	
    The system allows the user to retrieve and withdraw a specific
	record from the genetic database based on the ID.Ref from the
	metadata generated during the submission process.

C. DATABASE INITIALISATION, MODIFICATION, AND SAVING

The genetic database initialises using a pickle file to perform CRUD operations. The data for the record and the entire genetic database is stored in a pickle file, allowing users to create a new file, store information in the file, modify information, and save the data. Before starting the process, the record will be empty. The function will then create a new empty file to record submitted data by users. After each iteration of CRUD operation, any changes will be saved before exiting the program.

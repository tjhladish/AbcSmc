#ifndef SQDB_SQDB_H
#define SQDB_SQDB_H

#include <string>

#include "sqlite3.h"

#ifdef _WIN32
#  include <tchar.h>
#  define SQDB_MAKE_TEXT(x) _TEXT(x)
#  define SQDB_STRLEN _tcslen
#  define SQDB_STRDUP _tcsdup
#else
#  define SQDB_MAKE_TEXT(x) (x) 
#  define SQDB_STRLEN strlen
#  define SQDB_STRDUP strdup
#endif

#if !defined(SQDB_UTF16) && !defined(SQDB_UTF8)
#  ifdef _WIN32
#    if defined(UNICODE) || defined(_UNICODE)
#      define SQDB_UTF16
#    else
#      define SQDB_UTF8
#    endif
#  else
#    define SQDB_UTF8
#  endif
#endif

#ifdef SQDB_UTF8
#  define SQDB_CHAR char 
#  define SQDB_STD_STRING std::string
#endif

#ifdef SQDB_UTF16
#  define SQDB_CHAR TCHAR 
#  define SQDB_STD_STRING std::wstring
#endif

namespace sqdb
{

class Exception
{
public:
  Exception(sqlite3* db);

  Exception(sqlite3* db, int errorCode);

  Exception(const SQDB_CHAR* errorMsg);

  ~Exception();

  int GetErrorCode() const;

  const SQDB_CHAR* GetErrorMsg() const;
private:
  int m_errorCode;
  SQDB_CHAR* m_errorMsg;
};

#define CHECK(db, returnCode) \
  if ( (returnCode) != SQLITE_OK ) throw Exception(db, returnCode) 

class RefCount
{
protected:
  RefCount();

  RefCount(const RefCount& x);
  RefCount& operator=(const RefCount& x);

  void IncRef();
  unsigned DecRef();

private:
  unsigned* m_refCount;
};

class Blob : public RefCount
{
public:
  Blob(const void* data, int size);

  Blob(const Blob& x);
  Blob& operator=(const Blob& x);

  int GetSize() const;
  const char* GetData() const;

  ~Blob();

private:
  char* m_data;
  int m_size;
};

class Convertor
{
public:
  Convertor(sqlite3* db, sqlite3_stmt* stmt, int field);

  operator int() const;
  operator long long() const;
  operator double() const;
  operator SQDB_STD_STRING() const;
  operator const SQDB_CHAR*() const;
  operator Blob() const;

  int GetInt() const;
  long long GetLongLong() const;
  double GetDouble() const;
  SQDB_STD_STRING GetString() const;
  const SQDB_CHAR* GetText() const;
  Blob GetBlob() const;

private:
  sqlite3* m_db;
  sqlite3_stmt* m_stmt;
  int m_field;
};

class Statement : public RefCount
{
public:
  Statement(sqlite3* db, sqlite3_stmt* stmt);

  Statement(const Statement& x);
  Statement& operator=(const Statement& x);

  bool Next();
  Convertor GetField(int field) const;

  template<class T>
  void Bind(int i, const T& value)
  {
    if ( m_needReset ) 
      Reset();
    DoBind(i, value);
  }

  void BindBlob(int i, const void* value, int n);
  void BindNull(int i);

  ~Statement();

private:
  void DoBind(int i, int value); 
  void DoBind(int i, long long value); 
  void DoBind(int i, double value);
  void DoBind(int i, const SQDB_STD_STRING& value);
  void DoBind(int i, const SQDB_CHAR* value);

  // Bind blob.
  void DoBind(int i, const void* value, int n);

  // Bind null.
  void DoBind(int i);

  // Reset binders so that new values can be bound.
  void Reset();

  sqlite3* m_db;
  sqlite3_stmt* m_stmt;
  bool m_needReset;
};

class QueryStr
{
public:
  QueryStr();

  const SQDB_CHAR* Format(const SQDB_CHAR* fmt, ...);

  const SQDB_CHAR* Get() const;

  ~QueryStr();

private:
  SQDB_CHAR* m_buf;
};

class Db : public RefCount
{
public:
  Db(const SQDB_CHAR* fileName);

  void BeginTransaction();
  void CommitTransaction();
  void RollbackTransaction();

  bool TableExists(const SQDB_CHAR* tableName);
  Statement Query(const SQDB_CHAR* queryStr);
  long long LastId();

  Db(const Db& x);
  Db& operator=(const Db& x);

  ~Db();

private:
  sqlite3* m_db;
};

}

#endif


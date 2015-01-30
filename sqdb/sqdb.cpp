#include <cassert>
#include <cstdio>
#include <cstring>
#include <memory>

#include <stdlib.h>

#include "sqdb.h"

using namespace sqdb;

Exception::Exception(sqlite3* db)
{
  m_errorCode = sqlite3_errcode(db); 
  SQDB_CHAR* c = (SQDB_CHAR*)
#ifdef SQDB_UTF8
    sqlite3_errmsg
#else 
    sqlite3_errmsg16
#endif
    (db);
  m_errorMsg = SQDB_STRDUP(c);
}

Exception::Exception(sqlite3* db, int errorCode)
: m_errorCode(errorCode)
{
  SQDB_CHAR* c = (SQDB_CHAR*)
#ifdef SQDB_UTF8
    sqlite3_errmsg
#else 
    sqlite3_errmsg16
#endif
    (db);
  m_errorMsg = SQDB_STRDUP(c);
}

Exception::Exception(const SQDB_CHAR* errorMsg)
: m_errorCode(-1)
{
  m_errorMsg = SQDB_STRDUP(errorMsg);
}

Exception::~Exception()
{
  free(m_errorMsg);
}

int Exception::GetErrorCode() const
{
  return m_errorCode;
}

const SQDB_CHAR* Exception::GetErrorMsg() const
{
  return m_errorMsg;
}

RefCount::RefCount()
: m_refCount(NULL)
{
}

RefCount::RefCount(const RefCount& x)
: m_refCount(x.m_refCount)
{
}

RefCount& RefCount::operator=(const RefCount& x)
{
  if ( this != &x )
  {
    m_refCount = x.m_refCount;
  }
  return *this;
}

void RefCount::IncRef()
{
  if ( !m_refCount ) 
  {
    m_refCount = new unsigned;
    *m_refCount = 0;
  }
  ++*m_refCount;
}

unsigned RefCount::DecRef()
{
  assert(m_refCount);
  unsigned value = --*m_refCount;
  if ( value == 0 ) 
  {
    delete m_refCount;
    m_refCount = NULL;
  }
  return value;
}

Blob::Blob(const void* data, int size)
: m_size(size)
{
  m_data = new char[size];
  std::uninitialized_copy((char*)data, (char*)data + size, m_data); 
  IncRef();
}

Blob::Blob(const Blob& x)
: RefCount(x), m_data(x.m_data), m_size(x.m_size)
{
  IncRef();
}

Blob& Blob::operator=(const Blob& x)
{
  if ( this != &x )
  {
    RefCount::operator=(x);
    IncRef();
    m_data = x.m_data;
    m_size = x.m_size;
  }
  return *this;
}

int Blob::GetSize() const
{
  return m_size;
}

const char* Blob::GetData() const
{
  return m_data;
}

Blob::~Blob()
{
  if ( DecRef() == 0 ) 
  {
    delete[] m_data;
  }
}

Convertor::Convertor(sqlite3* db, sqlite3_stmt* stmt, int field)
: m_db(db), m_stmt(stmt), m_field(field)
{
}

Convertor::operator int() const
{
  return GetInt();
}

Convertor::operator long long() const
{
  return GetLongLong();
}

Convertor::operator double() const
{
  return GetDouble();
}

Convertor::operator SQDB_STD_STRING() const
{
  return GetString();
}

Convertor::operator const SQDB_CHAR*() const
{
  return GetText();
}

Convertor::operator Blob() const
{
  return GetBlob();
}

int Convertor::GetInt() const
{
  assert(m_stmt);
  return sqlite3_column_int(m_stmt, m_field);
}

long long Convertor::GetLongLong() const
{
  assert(m_stmt);
  return sqlite3_column_int64(m_stmt, m_field);
}

double Convertor::GetDouble() const
{
  assert(m_stmt);
  return sqlite3_column_double(m_stmt, m_field);
}

SQDB_STD_STRING Convertor::GetString() const
{
  assert(m_stmt);
  const SQDB_CHAR* result = (const SQDB_CHAR*)
#ifdef SQDB_UTF8
  sqlite3_column_text
#else
  sqlite3_column_text16
#endif
  (m_stmt, m_field);
  return result;
}

const SQDB_CHAR* Convertor::GetText() const
{
  assert(m_stmt);
  const SQDB_CHAR* result = (const SQDB_CHAR*)
#ifdef SQDB_UTF8
  sqlite3_column_text
#else
  sqlite3_column_text16
#endif
  (m_stmt, m_field);
  return result;
}

Blob Convertor::GetBlob() const
{
  assert(m_stmt);
  const void* data = sqlite3_column_blob(m_stmt, m_field);
  int size = sqlite3_column_bytes(m_stmt, m_field);
  return Blob(data, size);
}

Statement::Statement(sqlite3* db, sqlite3_stmt* stmt)
: RefCount(), m_db(db), m_stmt(stmt), m_needReset(false)
{
  IncRef();
}

Statement::Statement(const Statement& x)
: RefCount(x), m_db(x.m_db), m_stmt(x.m_stmt), m_needReset(false)
{
  IncRef();
}

Statement& Statement::operator=(const Statement& x)
{
  if ( this != &x )
  {
    RefCount::operator=(x);
    IncRef();
    m_db = x.m_db;
    m_stmt = x.m_stmt;
    m_needReset = x.m_needReset;
  }
  return *this;
}

bool Statement::Next()
{
  assert(m_stmt);
  int ret = sqlite3_step(m_stmt);
  m_needReset = true;
  if ( ret == SQLITE_DONE )
  {
    return false;
  }
  else if ( ret == SQLITE_ROW ) 
  {
    return true;
  }
  else
  {
    CHECK(m_db, ret);
  }
}

Convertor Statement::GetField(int field) const
{
  return Convertor(m_db, m_stmt, field);
}

Statement::~Statement()
{
  if ( DecRef() == 0 ) 
  {
    sqlite3_finalize(m_stmt);
  }
}

void Statement::BindBlob(int i, const void* value, int n)
{
  if ( m_needReset ) 
    Reset();
  DoBind(i, value, n);
}

void Statement::BindNull(int i)
{ 
  if ( m_needReset ) 
    Reset();
  DoBind(i);
}

void Statement::DoBind(int i, int value)
{
  const int ret = sqlite3_bind_int(m_stmt, i, value);
  CHECK(m_db, ret);
}

void Statement::DoBind(int i, long long value)
{
  const int ret = sqlite3_bind_int64(m_stmt, i, value);
  CHECK(m_db, ret);
}

void Statement::DoBind(int i, double value)
{
  const int ret = sqlite3_bind_double(m_stmt, i, value);
  CHECK(m_db, ret);
}

void Statement::DoBind(int i, const SQDB_STD_STRING& value)
{
  const int ret = 
#ifdef SQDB_UTF8
  sqlite3_bind_text
#else
  sqlite3_bind_text16
#endif
  (m_stmt, i, value.c_str(), value.size() * sizeof(SQDB_STD_STRING::value_type), SQLITE_TRANSIENT);
  CHECK(m_db, ret);
}

void Statement::DoBind(int i, const SQDB_CHAR* value)
{
  const int len = SQDB_STRLEN(value);

  const int ret = 
#ifdef SQDB_UTF8
  sqlite3_bind_text
#else
  sqlite3_bind_text16
#endif
  (m_stmt, i, value, len * sizeof(SQDB_CHAR), SQLITE_TRANSIENT);

  CHECK(m_db, ret);
}

void Statement::DoBind(int i, const void* value, int n)
{
  const int ret = sqlite3_bind_blob(m_stmt, i, value, n, SQLITE_TRANSIENT);
  CHECK(m_db, ret);
}

void Statement::DoBind(int i)
{
  const int ret = sqlite3_bind_null(m_stmt, i);
  CHECK(m_db, ret);
}

void Statement::Reset()
{
  assert(m_needReset);
  assert(m_stmt);

  sqlite3_reset(m_stmt);
  m_needReset = false;
}

QueryStr::QueryStr()
: m_buf(NULL)
{
}

const SQDB_CHAR* QueryStr::Format(const SQDB_CHAR* fmt, ...)
{
  va_list va;
  va_start(va, fmt);
#ifdef SQDB_UTF8
  sqlite3_free(m_buf);
  m_buf = sqlite3_vmprintf(fmt, va);
#else
  free(m_buf);
  int len = _vscwprintf(fmt, va) + 1;
  m_buf = (SQDB_CHAR*)malloc(len * sizeof(SQDB_CHAR));
  vswprintf(m_buf, len, fmt, va);
#endif

  va_end(va);

  return m_buf;
}

const SQDB_CHAR* QueryStr::Get() const
{
  return m_buf;
}

QueryStr::~QueryStr()
{
#ifdef SQDB_UTF8
  sqlite3_free(m_buf);
#else
  free(m_buf);
#endif
}

Db::Db(const SQDB_CHAR* fileName)
: RefCount()
{
#ifdef SQDB_UTF8
  const int ret = sqlite3_open(fileName, &m_db);
#else
  const int ret = sqlite3_open16(fileName, &m_db);
#endif
  CHECK(m_db, ret);

  IncRef();
}

void Db::BeginTransaction()
{
  Query(SQDB_MAKE_TEXT("BEGIN;")).Next();
}

void Db::CommitTransaction()
{
  Query(SQDB_MAKE_TEXT("COMMIT;")).Next();
}

void Db::RollbackTransaction()
{
  Query(SQDB_MAKE_TEXT("ROLLBACK;")).Next();
}

bool Db::TableExists(const SQDB_CHAR* tableName)
{
  QueryStr str;
  Statement s = 
    Query(
      str.Format(SQDB_MAKE_TEXT("select count(*) from sqlite_master where type='table' and name='%s';"), tableName));
  s.Next();
  const int count = s.GetField(0);
  return count > 0;
}

Statement Db::Query(const SQDB_CHAR* queryStr)
{
  sqlite3_stmt* stmt = NULL;
#ifdef SQDB_UTF8
  const int ret = sqlite3_prepare(m_db, queryStr, -1, &stmt, NULL);
#else
  const int ret = sqlite3_prepare16(m_db, queryStr, -1, &stmt, NULL);
#endif
  CHECK(m_db, ret);

  return Statement(m_db, stmt);
}

long long Db::LastId()
{
  long long ret = sqlite3_last_insert_rowid(m_db);
  return ret;
}

Db::Db(const Db& x)
: RefCount(x),
  m_db(x.m_db)
{
  IncRef();
}

Db& Db::operator=(const Db& x)
{
  if ( this != &x )
  {
    RefCount::operator=(x);

    IncRef();
    m_db = x.m_db;
  }
  return *this;
}

Db::~Db()
{
  if ( DecRef() == 0 ) 
  {
    sqlite3_close(m_db);
  }
}


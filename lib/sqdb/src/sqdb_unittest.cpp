#include <gtest/gtest.h>
#include "sqdb.h"

using namespace sqdb;

class SqDbTest : public ::testing::Test
{
protected:
  virtual void SetUp()
  {
    m_db = new Db("test.sdb");
    m_db->Query("create table t1(id int, str varchar(200))").Next();
    m_db->Query("create table t2(b blob)").Next();
  }

  virtual void TearDown()
  {
    delete m_db;
    unlink("test.sdb");
  }

  Db* m_db;
};

TEST_F(SqDbTest, TableExistsTest)
{
  EXPECT_TRUE(m_db->TableExists("t1"));
  EXPECT_TRUE(m_db->TableExists("t2"));
  EXPECT_FALSE(m_db->TableExists("hlajsflajsf"));
}

TEST_F(SqDbTest, InsertionTest)
{
  m_db->Query("insert into t1 values(1, 'jljjj')").Next();
  m_db->Query("insert into t1 values(2, 'jhhhh')").Next();
  Statement s = m_db->Query("select count(*) from t1"); 
  ASSERT_TRUE(s.Next());
  EXPECT_EQ(2, (int)s.GetField(0));
  ASSERT_FALSE(s.Next());
}

TEST_F(SqDbTest, ResetBindTest)
{
  Statement i = m_db->Query("insert into t1 values(?, ?)");
  i.Bind(1, 1);
  i.Bind(2, "ssss");
  i.Next();
  i.Bind(1, 2);
  i.Bind(2, "ljlkj");
  i.Next();
}

TEST_F(SqDbTest, BindInsertionTest)
{
  Statement i = m_db->Query("insert into t1 values(?, ?)");
  i.Bind(1, 1);
  i.Bind(2, "ssss");
  i.Next();
  i.Bind(1, 2);
  i.Bind(2, "ljlkj");
  i.Next();
  Statement s = m_db->Query("select count(*) from t1"); 
  ASSERT_TRUE(s.Next());
  EXPECT_EQ(2, (int)s.GetField(0));
  ASSERT_FALSE(s.Next());
}

TEST_F(SqDbTest, BlobTest)
{
  const size_t s = 1024;
  int* a = (int*)malloc(s * sizeof(int)); 
  for ( int i = 0; i < s; ++i, a[i] = i );
  Statement i = m_db->Query("insert into t2 (b) values(?)");
  i.BindBlob(1, a, s * sizeof(int));
  i.Next();
  Statement st = m_db->Query("select count(*) from t2");
  st.Next();
  EXPECT_EQ(1, (int)st.GetField(0));
  free(a);
}

TEST_F(SqDbTest, BlobTest2)
{
  const size_t s = 1024;
  int* a = (int*)malloc(s * sizeof(int)); 
  for ( int i = 0; i < s; a[i] = i, ++i );
  Statement i = m_db->Query("insert into t2 (b) values(?)");
  i.BindBlob(1, a, s * sizeof(int));
  i.Next();
  Statement st = m_db->Query("select b from t2");
  st.Next();
  Blob b = st.GetField(0);
  const int* d = (const int*)b.GetData();
  for ( int i = 0; i < s; ++i )
  {
    EXPECT_EQ(i, d[i]);
  }
  EXPECT_EQ(s * sizeof(int), b.GetSize());
  free(a);
}

TEST_F(SqDbTest, InsertionTransactionTest)
{
  try
  {
    m_db->BeginTransaction();
    // Wrong query.
    m_db->Query("insert into t1 values(1, 2, 'jljjj')").Next();
    m_db->CommitTransaction();
  }
  catch ( const sqdb::Exception& e )
  {
    m_db->RollbackTransaction();
  }
  Statement s = m_db->Query("select count(*) from t1"); 
  ASSERT_TRUE(s.Next());
  EXPECT_EQ(0, (int)s.GetField(0));
  ASSERT_FALSE(s.Next());
}

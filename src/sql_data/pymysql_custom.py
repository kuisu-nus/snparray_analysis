import mysql.connector
from mysql.connector import Error

def create_database():
    """创建数据库"""
    try:
        # 先连接到MySQL（不指定数据库）
        connection = mysql.connector.connect(
            host='localhost',
            user='root',
            password='sukui1016',
            auth_plugin='mysql_native_password'
        )
        
        if connection.is_connected():
            cursor = connection.cursor()
            
            # 创建数据库
            cursor.execute("CREATE DATABASE IF NOT EXISTS student_management")
            print("数据库 'student_management' 创建成功")
            
            cursor.close()
            connection.close()
            
    except Error as e:
        print(f"创建数据库错误: {e}")

def create_connection():
    """创建数据库连接"""
    try:
        connection = mysql.connector.connect(
            host='localhost',
            database='student_management',  # 现在这个数据库应该存在了
            user='root',
            password='sukui1016',
            auth_plugin='mysql_native_password'
        )
        if connection.is_connected():
            print("成功连接到MySQL数据库")
            return connection
    except Error as e:
        print(f"连接错误: {e}")
        return None

def create_tables(connection):
    """创建学生和课程表"""
    create_students_table = """
    CREATE TABLE IF NOT EXISTS students (
        student_id INT AUTO_INCREMENT PRIMARY KEY,
        name VARCHAR(100) NOT NULL,
        email VARCHAR(100) UNIQUE NOT NULL,
        age INT,
        major VARCHAR(100),
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
    """
    
    create_courses_table = """
    CREATE TABLE IF NOT EXISTS courses (
        course_id INT AUTO_INCREMENT PRIMARY KEY,
        course_name VARCHAR(100) NOT NULL,
        instructor VARCHAR(100),
        credits INT DEFAULT 3
    )
    """
    
    try:
        cursor = connection.cursor()
        cursor.execute(create_students_table)
        cursor.execute(create_courses_table)
        connection.commit()
        print("表创建成功")
    except Error as e:
        print(f"创建表错误: {e}")
    finally:
        cursor.close()


def insert_student(connection, name, email, age, major):
    """插入学生数据"""
    query = """
    INSERT INTO students (name, email, age, major)
    VALUES (%s, %s, %s, %s)
    """
    
    try:
        cursor = connection.cursor()
        cursor.execute(query, (name, email, age, major))
        connection.commit()
        print(f"学生 {name} 添加成功，ID: {cursor.lastrowid}")
        return cursor.lastrowid
    except Error as e:
        print(f"插入错误: {e}")
        connection.rollback()
        return None
    finally:
        cursor.close()

def insert_course(connection, course_name, instructor, credits=3):
    """插入课程数据"""
    query = """
    INSERT INTO courses (course_name, instructor, credits)
    VALUES (%s, %s, %s)
    """
    
    try:
        cursor = connection.cursor()
        cursor.execute(query, (course_name, instructor, credits))
        connection.commit()
        print(f"课程 {course_name} 添加成功")
        return cursor.lastrowid
    except Error as e:
        print(f"插入错误: {e}")
        connection.rollback()
        return None
    finally:
        cursor.close()


def get_all_students(connection):
    """获取所有学生"""
    query = "SELECT * FROM students"
    
    try:
        cursor = connection.cursor(dictionary=True)  # 返回字典格式
        cursor.execute(query)
        students = cursor.fetchall()
        return students
    except Error as e:
        print(f"查询错误: {e}")
        return []
    finally:
        cursor.close()

def get_student_by_id(connection, student_id):
    """根据ID获取学生"""
    query = "SELECT * FROM students WHERE student_id = %s"
    
    try:
        cursor = connection.cursor(dictionary=True)
        cursor.execute(query, (student_id,))
        student = cursor.fetchone()
        return student
    except Error as e:
        print(f"查询错误: {e}")
        return None
    finally:
        cursor.close()

def search_students(connection, name_pattern=None, major=None):
    """条件查询学生"""
    query = "SELECT * FROM students WHERE 1=1"
    params = []
    
    if name_pattern:
        query += " AND name LIKE %s"
        params.append(f"%{name_pattern}%")
    
    if major:
        query += " AND major = %s"
        params.append(major)
    
    try:
        cursor = connection.cursor(dictionary=True)
        cursor.execute(query, params)
        students = cursor.fetchall()
        return students
    except Error as e:
        print(f"查询错误: {e}")
        return []
    finally:
        cursor.close()


def update_student(connection, student_id, name=None, email=None, age=None, major=None):
    """更新学生信息"""
    # 动态构建更新语句
    update_fields = []
    params = []
    
    if name:
        update_fields.append("name = %s")
        params.append(name)
    if email:
        update_fields.append("email = %s")
        params.append(email)
    if age is not None:
        update_fields.append("age = %s")
        params.append(age)
    if major:
        update_fields.append("major = %s")
        params.append(major)
    
    if not update_fields:
        print("没有提供更新字段")
        return False
    
    params.append(student_id)
    query = f"UPDATE students SET {', '.join(update_fields)} WHERE student_id = %s"
    
    try:
        cursor = connection.cursor()
        cursor.execute(query, params)
        connection.commit()
        if cursor.rowcount > 0:
            print(f"学生 ID {student_id} 更新成功")
            return True
        else:
            print("没有找到要更新的学生")
            return False
    except Error as e:
        print(f"更新错误: {e}")
        connection.rollback()
        return False
    finally:
        cursor.close()

def delete_student(connection, student_id):
    """删除学生"""
    query = "DELETE FROM students WHERE student_id = %s"
    
    try:
        cursor = connection.cursor()
        cursor.execute(query, (student_id,))
        connection.commit()
        if cursor.rowcount > 0:
            print(f"学生 ID {student_id} 删除成功")
            return True
        else:
            print("没有找到要删除的学生")
            return False
    except Error as e:
        print(f"删除错误: {e}")
        connection.rollback()
        return False
    finally:
        cursor.close()

def transactional_operation(connection):
    """演示事务操作"""
    try:
        cursor = connection.cursor()
        
        # 开始事务
        connection.start_transaction()
        
        # 执行多个操作
        cursor.execute("INSERT INTO students (name, email, age, major) VALUES (%s, %s, %s, %s)",
                      ("事务测试", "transaction@test.com", 22, "Computer Science"))
        
        cursor.execute("UPDATE students SET age = age + 1 WHERE name = %s", ("事务测试",))
        
        # 提交事务
        connection.commit()
        print("事务执行成功")
        
    except Error as e:
        # 回滚事务
        connection.rollback()
        print(f"事务执行失败，已回滚: {e}")
    finally:
        cursor.close()


class DatabaseManager:
    """数据库管理类，封装常用操作和错误处理"""
    
    def __init__(self):
        self.connection = None
    
    def __enter__(self):
        self.connection = create_connection()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.connection and self.connection.is_connected():
            self.connection.close()
            print("数据库连接已关闭")
    
    def execute_query(self, query, params=None, fetch=False):
        """执行查询并处理错误"""
        try:
            cursor = self.connection.cursor(dictionary=True)
            cursor.execute(query, params or ())
            
            if fetch:
                result = cursor.fetchall()
            else:
                result = cursor.lastrowid
            
            self.connection.commit()
            return result
            
        except Error as e:
            print(f"数据库错误: {e}")
            self.connection.rollback()
            return None
        finally:
            if 'cursor' in locals():
                cursor.close()
    
    def batch_insert_students(self, students_data):
        """批量插入学生数据"""
        query = """
        INSERT INTO students (name, email, age, major)
        VALUES (%s, %s, %s, %s)
        """
        
        try:
            cursor = self.connection.cursor()
            cursor.executemany(query, students_data)
            self.connection.commit()
            print(f"成功插入 {cursor.rowcount} 条记录")
            return cursor.rowcount
        except Error as e:
            print(f"批量插入错误: {e}")
            self.connection.rollback()
            return 0
        finally:
            cursor.close()

def main():
    """主函数演示完整流程"""
    
    # 使用上下文管理器确保连接正确关闭
    with DatabaseManager() as db:
        if not db.connection:
            print("无法连接到数据库")
            return
        
        # 创建表
        create_tables(db.connection)
        
        # 插入单个学生
        student_id = insert_student(db.connection, "张三", "zhangsan@email.com", 20, "计算机科学")
        
        # 批量插入学生
        students_data = [
            ("李四", "lisi@email.com", 21, "数学"),
            ("王五", "wangwu@email.com", 22, "物理"),
            ("赵六", "zhaoliu@email.com", 19, "计算机科学")
        ]
        db.batch_insert_students(students_data)
        
        # 插入课程
        insert_course(db.connection, "Python编程", "张老师", 4)
        insert_course(db.connection, "数据库原理", "李老师", 3)
        
        # 查询所有学生
        print("\n所有学生:")
        students = get_all_students(db.connection)
        for student in students:
            print(f"ID: {student['student_id']}, 姓名: {student['name']}, 专业: {student['major']}")
        
        # 条件查询
        print("\n计算机科学专业的学生:")
        cs_students = search_students(db.connection, major="计算机科学")
        for student in cs_students:
            print(f"姓名: {student['name']}, 邮箱: {student['email']}")
        
        # 更新学生信息
        if student_id:
            update_student(db.connection, student_id, age=21, major="软件工程")
        
        # 演示事务
        transactional_operation(db.connection)
        
        # 最终查询结果
        print("\n最终学生列表:")
        final_students = get_all_students(db.connection)
        for student in final_students:
            print(f"ID: {student['student_id']}, 姓名: {student['name']}, 年龄: {student['age']}, 专业: {student['major']}")

if __name__ == "__main__":
    main()

if __name__ == '__main__':
    # 先创建数据库，再连接
    main()
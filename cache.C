#include <iostream>
#include <functional>

class CVU{
 public: 
  double foo() const {return 1.;}
};

class Event {
 public: 
  Event() : m_cvu(CVU()) {};
  CVU m_cvu;
  double foo() const {return 0.;}
};

class AbsVar {
 public:
  //typedef std::function<double(const T&)> PointerToTFunction;
  //PointerToTFunction m_pointer_to_GetValue;
  //AbsVar(PointerToTFunction p) : m_pointer_to_GetValue(p) {}
  virtual double GetValue (const Event& e) const = 0;
};

class EventVar : public AbsVar{
 public: 
  typedef std::function<double(const Event&)> PointerToTFunction;
  PointerToTFunction m_pointer_to_GetValue;
  EventVar(PointerToTFunction p) : m_pointer_to_GetValue(p) {};
  double GetValue (const Event& t) const { 
    return m_pointer_to_GetValue(t); 
  }
};

class CVUVar : public AbsVar {
 public: 
  typedef std::function<double(const CVU&)> PointerToTFunction;
  PointerToTFunction m_pointer_to_GetValue;
  CVUVar(PointerToTFunction p) : m_pointer_to_GetValue(p) {};
  double GetValue (const Event& t) const { 
    return m_pointer_to_GetValue(t.m_cvu); 
  }
};

int cache() {
  EventVar ev(&Event::foo);
  CVUVar cv(&CVU::foo);
  const Event e;
  const CVU c;
  std::vector<AbsVar*> vec = {&ev, &cv};
  for (auto i : vec) std::cout << i->GetValue(e) << "\n";
  return 0;
}
/*
  struct X{
      void show() {
          std::cout << "X\n";
      }
      //stuff
  }x1, x2, x3;

  struct Widget{
      void show() {
          std::cout << "Widget\n";
      }
      //stuff
  }w4, w5, w6;

  struct Toolbar{
      void show()
      {
          std::cout << "Toolbar\n";
      }
      //stuff
  }t1, t2, t3;


  template<class...Objects>
  void call_show(Objects&&...objects)
  {
      using expand = int[];
      (void) expand { 0, ((void)objects.show(), 0)... };
  }

  template<class Functor, class...Objects>
  void for_all(Functor&& f, Objects&&... objects)
  {
      using expand = int[];
      (void) expand { 0, (f(std::forward<Objects>(objects)), 0)... };

  }

  auto cache() -> int
  {
      call_show(x3, w4, w5, t1);
      for_all([](auto& thing) { thing.show(); }, x3, w4, w5, t1);
      return 0;
  }
*/

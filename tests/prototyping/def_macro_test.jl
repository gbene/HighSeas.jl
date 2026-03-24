macro def(name, definition)
  return quote
      macro $(esc(name))()
          esc($(Expr(:quote, definition)))
      end
  end
end

@def basic_person begin
    name::S
    age::Int
end

struct Person{S<:String}
    @basic_person
end


struct Citizen{S<:String}
    @basic_person
    nationality::String
end

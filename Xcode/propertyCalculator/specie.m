//
//  specie.m
//  propertyCalculator
//
//  Created by Niklas Nordin on 2011-01-05.
//  Copyright 2011 nequam. All rights reserved.
//

#import "specie.h"


@implementation specie

/*
-(void) init
{
	property_ = [[NSMutableArray alloc] init];
	function_ = [[NSMutableArray alloc] init];
}
*/
-(NSString *) name
{
	return name_;
}

-(void) setName:(NSString *)name
{
	name_ = name;
}
-(NSMutableArray *) properties
{
	return properties_;
}
@end
